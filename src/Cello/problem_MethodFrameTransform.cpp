// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFrameTransform.cpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-09-23
/// @brief    Implements the MethodFrameTransform class
///
/// The MethodFrameTransform method computes the average weighted velocity
/// (specified) component(s) where the weighting is performed using a specified
/// field. The field used for weighting should measure density, passive scalar
/// density or some other quantity measured per unit volume (the calculation
/// multiplies cell volume by the field values). It then modifies an attribute
/// of Block that tracks the current velocity of the reference frame (relative
/// to the frame when the velocity was initialized). It also updates the
/// "total_energy" field (since the kinetic energy changes) and the relevant
/// velocity fields.
///
/// This method also updates the origin_offset member of Block to indicate the
/// translation of the frame since the start of the simulation.

// Note that the comments and variable names all assume that the weight field
// is some kind of density. Thus they frequently refer to the weight field
// multiplied by cell volume as a "mass" and the product of this mass with
// velocity as a "momentum". We emphasize that as long as the weight field is
// simply a quantity per unit volume, this method can be used.

#include "problem.hpp"

// #define DEBUG_FRAME_TRANSFORM

#ifdef DEBUG_FRAME_TRANSFORM
#   define TRACE_FRAME_TRANSFORM CkPrintf ("%s:%d TRACE DEBUG_FRAME_TRANSFORM\n",__FILE__,__LINE__);
#else
#   define TRACE_FRAME_TRANSFORM /*   */
#endif

//----------------------------------------------------------------------

MethodFrameTransform::MethodFrameTransform
(bool component_transform[3], std::string weight_field,
 bool cycle_based_update, double update_start, double update_step,
 double weight_threshold, std::string threshold_type)
  : Method()
{
  cycle_based_update_ = cycle_based_update;
  prev_time_update_ = update_start - 0.5*update_step;

  ASSERT("MethodFrameTransform", "update_step must be >=0", update_step > 0.);
  cycle_start_ = (cycle_based_update_) ? (int) update_start :            0;
  cycle_step_  = (cycle_based_update_) ? (int)  update_step :            0;
  time_start_  = (cycle_based_update_) ?                 0. : update_start;
  time_step_   = (cycle_based_update_) ?                 0. :  update_step;

  // Copy component_transform entries
  int num_components = 0;
  for (int i=0; i<3; i++){
    component_transform_[i] = component_transform[i];
    num_components += (int)(component_transform_[i]);
  }

  // Check that the we will be performing transformations for at least velocity
  // component
  ASSERT("MethodFrameTransform",
	 ("Reference frame transformations are not being performed on any "
	  "velocity components"),
	 num_components > 0);

  // Check that the velocity components that we will be performing
  // transformations for actually make sense.
  int rank = cello::rank();
  if (rank == 1){
    ASSERT("MethodFrameTransform",
	   ("The y and z velocity components can't be transformed for a 1D "
	    "simulation."),
	   !(component_transform_[1] || component_transform_[2]));
  } else if (rank == 2){
    ASSERT("MethodFrameTransform",
	   ("The z velocity components can't be transformed for a 2D"
	    "simulation."), !component_transform_[2]);
  }

  FieldDescr * field_descr = cello::field_descr();

  // Check that velocity_x is a real field (important that it exists for
  // boundary conditions).
  ASSERT("MethodFrameTransform", "\"velocity_x\" field is not found.",
	 field_descr->is_field("velocity_x"));

  // Check that the weight_field is in fact a real, permanent field
  ASSERT1("MethodFrameTransform",
	  "\"%s\" is not the name of a permanent field",
	  field_descr->is_field(weight_field), weight_field.c_str());

  // Finally, store the weight_field value
  weight_field_ = weight_field;

  // we could probably reduce the ghost depth to zero (and possibly change the
  // synch type)
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
			     sync_id_method_frame_transform);
  refresh(ir)->add_field(field_descr->field_id(weight_field));
  if (field_descr->is_field("total_energy")){
    refresh(ir)->add_field(field_descr->field_id("total_energy"));
  }
  if (rank >= 1) refresh(ir)->add_field(field_descr->field_id("velocity_x"));
  if (rank >= 2) refresh(ir)->add_field(field_descr->field_id("velocity_y"));
  if (rank >= 3) refresh(ir)->add_field(field_descr->field_id("velocity_z"));

  weight_threshold_ = weight_threshold;

  // determine the type of threshold:
  std::string formatted(threshold_type.size(), ' ');
  std::transform(threshold_type.begin(), threshold_type.end(),
		 formatted.begin(), ::tolower);
  if ((formatted == "") || (formatted == "ignore")) {
    threshold_type_ = threshold_enum::ignore;
  } else if (formatted == "lower_limit") {
    threshold_type_ = threshold_enum::lower_limit;
  } else if (formatted == "upper_limit") {
    threshold_type_ = threshold_enum::upper_limit;
  } else {
    ERROR("MethodFrameTransform",
	  ("threshold_type must be an empty string, \"ignore\", "
	   "\"lower_limit\" or \"upper_limit\""));
  }
}

//----------------------------------------------------------------------

void MethodFrameTransform::pup(PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  PUParray(p,component_transform_,3);
  p|weight_field_;
  p|weight_threshold_;
  p|cycle_based_update_;
  p|cycle_start_;
  p|cycle_step_;
  p|time_start_;
  p|time_step_;
  p|prev_time_update_;
  // p|threshold_type_; results in errors on some systems (ex: using MPI)
  // therefore, the following is necessary
  if (p.isUnpacking()){
    int temp;
    p|temp;
    threshold_type_ =  static_cast<threshold_enum>(temp);
  } else {
    int temp = static_cast<int>(threshold_type_);
    p|temp;
  }
}

//----------------------------------------------------------------------

void MethodFrameTransform::compute( Block * block) throw()
{
  TRACE_FRAME_TRANSFORM;

  ASSERT("MethodFrameTransform::compute", "Not compatible with AMR",
	 block->level() == 0 && block->is_leaf());

  // Exit early if current time/cycle preceds first frame velocity update
  if (precede_first_transform_(block)){
    block->compute_done();
    return;
  }

  // update the current position of the block (from last cycle)
  double x_origin, y_origin, z_origin;
  block->origin_offset(&x_origin, &y_origin, &z_origin);
  double vx, vy, vz;
  block->frame_velocity(&vx,&vy,&vz);
  double dt = block->dt();
  block->set_origin_offset(x_origin + vx * dt,
			   y_origin + vy * dt,
			   z_origin + vz * dt);

  // check if frame velocity update is scheduled for current time/cycle. If an
  // update is scheduled, and scheduling is time-based, the following line also
  // records that the update will be made (needed for future schedule checks)
  bool scheduled_update = scheduled_velocity_update_(block);
  if (!scheduled_update) { // Exit early if a velocity update is not scheduled
    block->compute_done();
    return;
  }

  Field field = block->data()->field();

  // Reminder: we refer to the product of weight_field values and cell volumes
  // as "masses" and the product of these "masses" with velocity components as
  // "momentum" components

  // Compute the total mass AND momentum (only specified components) for the
  // weighted field on this block. The momentum components corresponding to
  // untracked velocity components are set to 0.
  double mass = 0;
  double momentum[3] = {0.,0.,0.};

  if (block->is_leaf()) {
    // Sum the compute mass and momentum from the local block
    precision_type precision = field_precision_(field);
    switch (precision) {
    case precision_single:
      { block_totals_<float>(block, mass, momentum); }
      break;
    case precision_double:
      { block_totals_<double>(block, mass, momentum); }
      break;
    case precision_extended80:
    case precision_extended96:
    case precision_quadruple:
      { block_totals_<long double>(block, mass, momentum); }
      break;
    }
  }

  // Assign the block an index based on it's location in the grid. This will
  // need to be modified once AMR and solvers are in use.
  int ix, iy, iz, nx, ny, nz;
  block->index_global(&ix, &iy, &iz, &nx, &ny, &nz);
  int index = ix + nx*(iy + ny*iz);

  // now, call the reduction to sum the momentum and mass over all blocks
  // (In principle, we could save bandwidth and cpu time by only performing the
  //  reduction on the relevant momentum components)

  double message_arr[5];
  message_arr[0] = (double)index;
  message_arr[1] = mass;
  for (int i=2; i<5; i++) { message_arr[i] = momentum[i-2]; }

  // Originally, we used did not include index in message_arr and used the
  // CkReduction::sum_double reduction type. Unfortunately, the order of
  // addition is undefined (depends on order that messages are received). To
  // make the results, we now collect all of the messages and sort them by
  // location before summing them.

  CkCallback cb(CkIndex_Block::p_method_frame_transform_end(NULL),
		block->proxy_array());

  block->contribute(5*sizeof(double), message_arr, CkReduction::concat, cb);
}

//----------------------------------------------------------------------

template <class T>
void MethodFrameTransform::block_totals_(Block * block,
					 double &mass,
					 double momentum[3])
  const throw()
{

  // Reminder: we refer to the product of weight_field values and cell volumes
  // as "masses" and the product of these "masses" with velocity components as
  // "momentum" components

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

  T* weight_vals = (T *) field.values(weight_field_);
  T * v3[3] = { (T*) (              field.values("velocity_x")),
		(T*) ((rank >= 2) ? field.values("velocity_y") : NULL),
		(T*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  // since volume within a block is constant, multiply by volume at the end
  mass = 0.;
  for (int j=0; j<3; j++) {momentum[j] = 0.;}


  // define lambda function to apply threshold.
  const double weight_threshold = weight_threshold_;
  const threshold_enum threshold_type = threshold_type_;

  auto satisfies_thresh = [=](double value)->bool{
    switch(threshold_type) {
    case(threshold_enum::ignore)      : { return true; }
    case(threshold_enum::lower_limit) : { return value >= weight_threshold; }
    case(threshold_enum::upper_limit) : { return value <= weight_threshold; }
    }
    return false;
  };


  // only contain material in active zone
  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {

	int i = ix + ndx*(iy + ndy*(iz));

	double density = (double) weight_vals[i];

	if (!satisfies_thresh(density)) { continue; }
	// add current density to running sum of densities
	mass += density;

	for (int j=0; j<3; j++){
	  // skip velocity components that aren't being transformed (this
	  // accounts for the rank of the simulation implicitly)
	  if (!component_transform_[j]) { continue; }
	  momentum[j] += (double) (v3[j][i]*density);
	}
      }
    }
  }

  // multiply mass and momentum by the volume
  double hx,hy,hz;
  block->data()->field_cell_width(&hx,&hy,&hz);
  double volume = hx*hy*hz;

  mass *= volume;
  for (int j=0; j<3; j++){momentum[j] *= volume;}

}
//----------------------------------------------------------------------

void Block::p_method_frame_transform_end(CkReductionMsg * msg)
{
  // This method is the entry method called after the CkReduction::sum_double,
  // specified in MethodFrameTransform::compute, is complete
  TRACE_FRAME_TRANSFORM;
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

template <class T>
void MethodFrameTransform::transform_field (Block * block,
					    T vx, T vy, T vz,
					    int ndx, int ndy, int ndz,
					    int nx,  int ny,  int nz,
					    int ix0, int iy0, int iz0,
					    int i_hist) throw()
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

template void MethodFrameTransform::transform_field
(Block * block,
 float vx, float vy, float vz,
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz,
 int ix0, int iy0, int iz0,
 int i_hist) throw();

template void MethodFrameTransform::transform_field
(Block * block,
 double vx, double vy, double vz,
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz,
 int ix0, int iy0, int iz0,
 int i_hist) throw();

template void MethodFrameTransform::transform_field
(Block * block,
 long double vx, long double vy, long double vz,
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz,
 int ix0, int iy0, int iz0,
 int i_hist) throw();

//----------------------------------------------------------------------

void MethodFrameTransform::compute_resume 
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_FRAME_TRANSFORM;

  // Reminder: we refer to the product of weight_field values and cell volumes
  // as "masses" and the product of these "masses" with velocity components as
  // "momentum" components

  // This method is called after the reduction is complete

  // To make our results deterministic, extract all the concatenated
  // reduction messages and sort them by the block index. We sort in place to
  // avoid memory allocation (the charm++ ampi implementation does this too)
  int n = msg->getSize()/(5*sizeof(double));
  double *values=(double *)msg->getData();

  for (int i = 0; i < n; i++){
    int index = (int) values[i*5]; // this was cast to double for convenience
    if (index != i){
      double temp_buf[5];
      for (int j = 0; j < 5; j++){ temp_buf[j] = values[i*5+j]; }

      while (index != i){
	index = (int) temp_buf[0];
	for (int j = 0; j < 5; j++){
	  double temp = values[index*5+j];
	  values[index*5+j] = temp_buf[j];
	  temp_buf[j] = temp;
	}
      }
    }
  }

  // Now sum up all mass and momentum
  double mass = 0.;
  double momentum[3] = {0., 0., 0.};
  for (int i=0; i<n; i++) {
    mass += values[i*5 + 1];
    for (int j=0; j<3; j++){
      momentum[j] += values[i*5 + j + 2];
    }
  }

  // compute the weighted velocity
  double frame_velocity[3];
  for (int i = 0; i<3; i++){
    if (mass > 0){
      frame_velocity[i] = momentum[i]/mass;
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

    // iterate over the data in the local block and compute the local total
    // "mass" (sum of weight fields satisfying threshold * cell_volume) and
    // local total "momentum"

    Field field = block->data()->field();

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
	MethodFrameTransform::transform_field (block, velocity[0],
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
	MethodFrameTransform::transform_field (block, velocity[0],
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
	MethodFrameTransform::transform_field (block, velocity[0],
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

precision_type MethodFrameTransform::field_precision_(Field &field)
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
	ASSERT("MethodFrameTransform::field_precision_",
	       "velocity and energy fields must have the same precision.",
	       out == p);
      }

    }
  }
  return out;
}

//----------------------------------------------------------------------

bool MethodFrameTransform::precede_first_transform_(Block * block)
  const throw()
{
  if (cycle_based_update_){
    return block->cycle() < cycle_start_;
  } else {
    return block->time()  > time_start_ ;
  }
}

//----------------------------------------------------------------------

bool MethodFrameTransform::scheduled_velocity_update_(Block * block)
  throw()
{
  if (cycle_based_update_){
    return ( ( (block->cycle() - cycle_start_) % cycle_step_ ) == 0 );
  } else {
    double cur_time = block->time();
    if (cur_time == prev_time_update_){ return true; }
    int cur_time_interval =
      (int) floor((cur_time - time_start_) / time_step_);
    int prev_time_interval =
      (int) floor((prev_time_update_ - time_start_) / time_step_);
    // a velocity update is only scheduled if current time interval is greater
    // than the time interval from the prior update during last update
    bool out = cur_time_interval > prev_time_interval;
    if (out){ prev_time_update_ = cur_time; }
    return out;
  }
}
