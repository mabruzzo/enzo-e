// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodScalarFrameTransform.cpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-09-23
/// @brief    Implements the MethodScalarFrameTransform class
///
/// The MethodScalarFrameTransform method computes the average
/// mass-weighted velocity component(s) (that can be freely specified) for a
/// specified passively advected scalar. It then modifies an attribute of Block
/// that tracks the current velocity of the reference frame (relative to the
/// frame when the velocity was initialized). It also updates the
/// "total_energy" field (since the kinetic energy changes) and the relevant
/// velocity fields
///
/// This requires that we perform a reduction

#include "problem.hpp"

// #define DEBUG_FRAME_TRANSFORM

#ifdef DEBUG_FRAME_TRANSFORM
#   define TRACE_FRAME_TRANSFORM CkPrintf ("%s:%d TRACE DEBUG_FRAME_TRANSFORM\n",__FILE__,__LINE__);
#else
#   define TRACE_FRAME_TRANSFORM /*   */
#endif

//----------------------------------------------------------------------

MethodScalarFrameTransform::MethodScalarFrameTransform
(bool component_transform[3], std::string passive_scalar,
 bool ignore_neg_scalar, int initial_cycle, int update_stride)
  : Method()
{
  ignore_neg_scalar_ = ignore_neg_scalar;
  initial_cycle_ = initial_cycle;
  ASSERT("MethodScalarFrameTransform", "update_stride must be >=0",
	 update_stride > 0);
  update_stride_ = update_stride;

  // Copy component_transform entries
  int num_components = 0;
  for (int i=0; i<3; i++){
    component_transform_[i] = component_transform[i];
    num_components += (int)(component_transform_[i]);
  }

  // Check that the we will be performing transformations for at least velocity
  // component
  ASSERT("MethodScalarFrameTransform",
	 ("Reference frame transformations are not being performed on any "
	  "velocity components"),
	 num_components > 0);

  // Check that the velocity components that we will be performing
  // transformations for actually make sense.
  int rank = cello::rank();
  if (rank == 1){
    ASSERT("MethodScalarFrameTransform",
	   ("The y and z velocity components can't be transformed for a 1D "
	    "simulation."),
	   !(component_transform_[1] || component_transform_[2]));
  } else if (rank == 2){
    ASSERT("MethodScalarFrameTransform",
	   ("The z velocity components can't be transformed for a 2D"
	    "simulation."), !component_transform_[2]);
  }

  FieldDescr * field_descr = cello::field_descr();

  // Check that velocity_x is a real field (important that it exists for
  // boundary conditions).
  ASSERT("MethodScalarFrameTransform", "\"velocity_x\" field is not found.",
	 field_descr->is_field("velocity_x"));

  // Check that the specified passively advected scalar is a real field
  ASSERT1("MethodScalarFrameTransform",
	  "\"%s\" is not the name of a permanent field",
	  field_descr->is_field(passive_scalar), passive_scalar.c_str());

  // Finally, store the passive_scalar value
  passive_scalar_ = passive_scalar;

  // we could probably reduce the ghost depth to zero (and possibly change the
  // synch type)
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     sync_id_method_scalar_frame_transform);
  refresh(ir)->add_field(field_descr->field_id(passive_scalar));
  if (field_descr->is_field("total_energy")){
    refresh(ir)->add_field(field_descr->field_id("total_energy"));
  }
  if (rank >= 1) refresh(ir)->add_field(field_descr->field_id("velocity_x"));
  if (rank >= 2) refresh(ir)->add_field(field_descr->field_id("velocity_y"));
  if (rank >= 3) refresh(ir)->add_field(field_descr->field_id("velocity_z"));
}

//----------------------------------------------------------------------

void MethodScalarFrameTransform::pup(PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  PUParray(p,component_transform_,3);
  p|passive_scalar_;
  p|ignore_neg_scalar_;
  p|initial_cycle_;
  p|update_stride_;
}

//----------------------------------------------------------------------

void MethodScalarFrameTransform::compute( Block * block) throw()
{
  TRACE_FRAME_TRANSFORM;

  int cur_cycle = block->cycle();

  // Exit Early if the frame velocity has not yet been modified
  if (cur_cycle < initial_cycle_){
    block->compute_done();
    return;
  }

  // update the current position of the block
  double x_origin, y_origin, z_origin;
  block->origin_offset(&x_origin, &y_origin, &z_origin);
  double vx, vy, vz;
  block->frame_velocity(&vx,&vy,&vz);
  double dt = block->dt();
  block->set_origin_offset(x_origin + vx * dt,
			   y_origin + vy * dt,
			   z_origin + vz * dt);

  // Exit Early if not scheduled to modify frame velocity
  if ( ((cur_cycle - initial_cycle_) % update_stride_) != 0 ) {
    block->compute_done();
    return;
  }


  Field field = block->data()->field();

  // Compute the total mass AND momentum (only specified components) for the
  // passive scalar on this block. The momentum components corresponding to
  // untracked velocity components are set to 0.
  double scalar_mass = 0;
  double scalar_momentum[3] = {0.,0.,0.};

  if (block->is_leaf()) {
    // Sum the compute scalar_mass and scalar_momentum from the local block
    precision_type precision = field_precision_(field);
    switch (precision) {
    case precision_single:
      { block_totals_<float>(block, scalar_mass, scalar_momentum); }
      break;
    case precision_double:
      { block_totals_<double>(block, scalar_mass, scalar_momentum); }
      break;
    case precision_extended80:
    case precision_extended96:
    case precision_quadruple:
      { block_totals_<long double>(block, scalar_mass, scalar_momentum); }
      break;
    }
  }

  // now, call the reduction to sum the momentum and mass over all blocks
  // (In principle, we could save bandwidth and cpu time by only performing the
  //  reduction on the relevant momentum components)

  double message_arr[4];
  message_arr[0] = scalar_mass;           message_arr[1] = scalar_momentum[0];
  message_arr[2] = scalar_momentum[1];    message_arr[3] = scalar_momentum[2];

  CkCallback cb(CkIndex_Block::p_method_scalar_frame_transform_end(NULL),
		block->proxy_array());

  block->contribute(4*sizeof(double),message_arr,
		    CkReduction::sum_double, cb);
}

//----------------------------------------------------------------------

template <class T>
void MethodScalarFrameTransform::block_totals_(Block * block,
					       double &scalar_mass,
					       double scalar_momentum[3])
  const throw()
{

  const int rank = cello::rank();

  Field field = block->data()->field();
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  T* scalar_d = (T *) field.values(passive_scalar_);
  T * v3[3] = { (T*) (              field.values("velocity_x")),
		(T*) ((rank >= 2) ? field.values("velocity_y") : NULL),
		(T*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  // since volume within a block is constant, multiply by volume at the end
  scalar_mass = 0.;
  for (int j=0; j<3; j++) {scalar_momentum[j] = 0.;}

  // only contain material in active zone
  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {

	int i = ix + ndx*(iy + ndy*(iz));

	T density = scalar_d[i];

	if (ignore_neg_scalar_ && density < 0) { continue; }
	// add current scalar density to running sum of scalar densities
	scalar_mass += (double) density;

	for (int j=0; j<3; j++){
	  // skip velocity components that aren't being transformed (this
	  // accounts for the rank of the simulation implicitly)
	  if (!component_transform_[j]) { continue; }
	  scalar_momentum[j] += (double) (v3[j][i]*density);
	}
      }
    }
  }

  // multiply scalar_mass and scalar_momentum by the volume
  double hx,hy,hz;
  block->data()->field_cell_width(&hx,&hy,&hz);
  double volume = hx*hy*hz;

  scalar_mass *= volume;
  for (int j=0; j<3; j++){scalar_momentum[j] *= volume;}

}
//----------------------------------------------------------------------

void Block::p_method_scalar_frame_transform_end(CkReductionMsg * msg)
{
  // This method is the entry method called after the CkReduction::sum_double,
  // specified in MethodScalarFrameTransform::compute, is complete
  TRACE_FRAME_TRANSFORM;
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void MethodScalarFrameTransform::transform_field (Block * block,
						  T vx, T vy, T vz,
						  int ndx, int ndy, int ndz,
						  int nx,  int ny,  int nz,
						  int ix0, int iy0, int iz0,
						  int i_hist)
  throw()
{
  // updates the relevant fields to the reference frame moving with
  // velocity = {vx,vy,vz} (with respect to the frame where fields were
  // initially measured)
  Field field = block->data()->field();

  const int rank = cello::rank();
  T * v3[3] =
    { (T*) (              field.values("velocity_x", i_hist)       ),
      (T*) ((rank >= 2) ? field.values("velocity_y", i_hist) : NULL),
      (T*) ((rank >= 3) ? field.values("velocity_z", i_hist) : NULL)  };

  T * te = NULL;
  if (field.is_field("total_energy")){
    te = (T*) field.values("total_energy", i_hist);
  }

  // identify velocity components not being transformed:
  bool skip_dim[3] = { (vx == 0.),
		       (vy == 0.) || (rank<2),
		       (vz == 0.) || (rank<3)};
  T frame_velocity[3] = {vx, vy, vz};

  // update the field values
  for (int ix=ix0; ix<ix0+nx; ix++) {
    for (int iy=iy0; iy<iy0+ny; iy++) {
      for (int iz=iz0; iz<iz0+nz; iz++) {
	int i = ix + ndx*(iy + ndy*(iz));

	T old_ke = 0; // ke from pre-transformed velocity components
	T new_ke = 0; // ke from transformed velocity components

	for (int j=0; j<3; j++){
	  // skip velocity components that aren't being transformed
	  if (skip_dim[j]){ continue; }

	  // before transforming velocity component, get old ke contribution
	  old_ke += 0.5*v3[j][i]*v3[j][i];

	  // transform velocity component
	  v3[j][i] -= frame_velocity[j];

	  // compute new ke contribution after transforming velocity component
	  new_ke += 0.5*v3[j][i]*v3[j][i];
	}

	if (te){
	  // update total energy
	  te[i] += new_ke - old_ke;
	}

      }
    }
  }
  
}


//----------------------------------------------------------------------

void MethodScalarFrameTransform::compute_resume 
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_FRAME_TRANSFORM;
  
  // This method is called after the reduction is complete
  double* result = (double *) msg->getData();
  double scalar_mass = result[0];
  double scalar_momentum[3] = {result[1], result[2], result[3]};

  // compute the average velocity of passively advected scalar
  double frame_velocity[3];
  for (int i = 0; i<3; i++){
    if (scalar_mass > 0){
      frame_velocity[i] = scalar_momentum[i]/scalar_mass;
    } else {
      frame_velocity[i] = 0.;
    }
    
  }

  // Update the velocity of the current frame measured with respect to the
  // frame when the gas was originally initialized (add the new frame velocity
  // to the old frame velocity)
  double vx,vy,vz;
  block->frame_velocity(&vx, &vy, &vz);
  block->set_frame_velocity(vx + frame_velocity[0],
			    vy + frame_velocity[1],
			    vz + frame_velocity[2]);

  // The rest of this only deals with leaf blocks:
  if (block->is_leaf()){

    Field field = block->data()->field();

    // iterate over the local grid and update the total_scalar_mass
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    int gx,gy,gz;
    field.ghost_depth (0,&gx,&gy,&gz);

    int ndx=nx+2*gx;
    int ndy=ny+2*gy;
    int ndz=nz+2*gz;

    // Note: Unclear if we should iterate over all history fields indices
    int i_hist = 0;

    // We are just iterating over the entire mesh. To just iterate over active
    // zone, don't modify nx,ny,nz and set ix0, iy0, iz0 to gx, gy, gz
    int ix0 = 0, iy0 = 0, iz0 = 0;

    nx = ndx;
    ny = ndy;
    nz = ndz;

    precision_type precision = field_precision_(field);

    switch (precision) {
    case precision_single:
      {
	float velocity[3];
	for (int j=0; j<3; j++) {velocity[j] =       (float)frame_velocity[j];}
	MethodScalarFrameTransform::transform_field (block, velocity[0],
						     velocity[1], velocity[2],
						     ndx, ndy, ndz,
						     nx,  ny,  nz,
						     ix0, iy0, iz0,
						     i_hist);
      }
      break;
    case precision_double:
      {
	double velocity[3];
	for (int j=0; j<3; j++) {velocity[j] =      (double)frame_velocity[j];}
	MethodScalarFrameTransform::transform_field (block, velocity[0],
						     velocity[1], velocity[2],
						     ndx, ndy, ndz,
						     nx,  ny,  nz,
						     ix0, iy0, iz0,
						     i_hist);
      }
      break;
    case precision_extended80:
    case precision_extended96:
    case precision_quadruple:
      {
	long double velocity[3];
	for (int j=0; j<3; j++) {velocity[j] = (long double)frame_velocity[j];}
	MethodScalarFrameTransform::transform_field (block, velocity[0],
						     velocity[1], velocity[2],
						     ndx, ndy, ndz,
						     nx,  ny,  nz,
						     ix0, iy0, iz0,
						     i_hist);
      }
      break;
    }
  }

  delete msg;
  block->compute_done();
}


//----------------------------------------------------------------------

precision_type MethodScalarFrameTransform::field_precision_(Field &field)
  const throw()
{
  bool identified_precision = false;
  precision_type out = precision_default; // initialize to satisfy compiler
  std::string names[4] = {"velocity_x", "velocity_y", "velocity_z",
			  "total_energy"};
  for (int i=0; i<4; i++){

    if (( (i <  3) && component_transform_[i]  ) ||
	( (i >= 3) && field.is_field(names[i]) )){
      precision_type p = field.precision(field.field_id(names[i]));
      if (!identified_precision){
	out = p;
	identified_precision = true;
      } else {
	ASSERT("MethodScalarFrameTransform::field_precision_",
	       "velocity and energy fields must have the same precision.",
	       out == p);
      }

    }
  }
  return out;
    
}
