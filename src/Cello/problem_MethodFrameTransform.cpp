// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFrameTransform.cpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-09-23
/// @brief    Implements the MethodFrameTransform class
///
/// Whenever MethodFrameTransform method is scheduled, it performs a reduction
/// on the all of the leaf blocks to determine a new frame velocity (measured
/// in the new reference frame) and updates the reference frame accordingly
///
/// The calculation of the new frame velocity is managed by the helper class
/// FrameTransformReductionMgr. See its declaration/definition for more details

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

#include "problem.hpp"

// #define DEBUG_FRAME_TRANSFORM

#ifdef DEBUG_FRAME_TRANSFORM
#   define TRACE_FRAME_TRANSFORM CkPrintf ("%s:%d TRACE DEBUG_FRAME_TRANSFORM\n",__FILE__,__LINE__);
#else
#   define TRACE_FRAME_TRANSFORM /*   */
#endif

template<typename T>
void select_val_(std::string arg_val, const std::map<std::string, T> &mapping,
                 T& out_val, std::string arg_name, std::string function_name){
  std::string formatted(arg_val.size(), ' ');
  std::transform(arg_val.begin(), arg_val.end(), formatted.begin(), ::tolower);
  auto it = mapping.find(formatted);
  if (it != mapping.cend()){
    out_val = it->second;
  } else {
    std::string error = arg_name + " must be ";
    std::size_t remaining = mapping.size() - 1;
    for (auto it = mapping.cbegin(); it != mapping.cend(); ++it){
      error += (it->first == "") ? "\"\"" : "\"" + it->first + "\"";
      if (remaining > 0){
        if (mapping.size() > 2) { error += ", "; };
        if (remaining == 1){ error += "or ";}
        remaining--;
      }
    }
    ERROR(function_name.c_str(), error.c_str());
  }
}

//----------------------------------------------------------------------

void check_target_downstream_dist_(frame_trans_reduce_enum reduction_type,
                                   const double (&target_downstream_dist)[3],
                                   const bool (&component_transform)[3])
{

  // as a sanity check, make sure that target_downstream_dist_ corresponds to
  // cell-centered locations on the coarsest grid
  const Config *config_ptr = cello::config();
  for (int i = 0; i<3; i++){
    if (reduction_type != frame_trans_reduce_enum::target_downstream_dist){
      ASSERT("FrameTransformReductionMgr",
             ("target_downstream_dist should not be specified for "
              "chosen reduction type."),
             target_downstream_dist[i] == 0);
    } else if ( (i >= cello::rank()) || (!component_transform[i]) ){
      ASSERT1("FrameTransformReductionMgr",
              ("target_downstream_dist[%d] should be 0 because no "
               "frame-tracking is performed along that axis."), i,
              target_downstream_dist[i] == 0);
    } else {

      double domain_width = (config_ptr->domain_upper[i] -
                             config_ptr->domain_lower[i]);
      int nroot_cells = config_ptr->mesh_root_size[i];

      ASSERT("FrameTransformReductionMgr",
             "target_downstream_dist must not be negative",
             target_downstream_dist[i] >= 0);
      ASSERT("FrameTransformReductionMgr",
             "target_downstream_dist_ must not exceed the domain width",
             target_downstream_dist[i] < domain_width);

      // if we use this, we need to be more careful about the size of
      // nroot_cells if it's close to 2^53, it will be problematic (the storage
      // of 2^53+1 in a double loses precision)
      double cell_width = domain_width/nroot_cells;
      double root_index = floor(target_downstream_dist[i] / cell_width);

      if ((cell_width*1e-6) < (cello::err_abs((root_index+0.5)*cell_width,
                                              target_downstream_dist[i]))){
        ERROR4("FrameTransformReductionMgr",
               ("Along axis-%d, the target downstream distance, %e, is not "
                "aligned with a cell-center. The nearest valid values are %e "
                " and %e"),
               i, target_downstream_dist[i],
               (root_index+0.5)*cell_width, (root_index+1.5)*cell_width);
      }
    }

  }
}

//----------------------------------------------------------------------

FrameTransformReductionMgr::FrameTransformReductionMgr
(std::string weight_field, double weight_threshold, std::string threshold_type,
 std::string reduction_type, const double (&target_downstream_dist)[3],
 const bool (&component_transform)[3]) throw()
{
  weight_field_ = weight_field;
  weight_threshold_ = weight_threshold;

  static const std::map<std::string, threshold_enum> thresh_map =
    {{"", threshold_enum::ignore},
     {"ignore", threshold_enum::ignore},
     {"lower_limit", threshold_enum::lower_limit},
     {"upper_limit", threshold_enum::upper_limit}};
  select_val_(threshold_type, thresh_map, threshold_type_, "threshold_type",
              "FrameTransformReductionMgr");

  static const std::map<std::string, frame_trans_reduce_enum> reduce_map =
    {{"weighted_average", frame_trans_reduce_enum::weighted_average},
     {"min", frame_trans_reduce_enum::min},
     {"min_zero_floor", frame_trans_reduce_enum::min_zero_floor},
     {"target_downstream_dist",
      frame_trans_reduce_enum::target_downstream_dist}};

  select_val_(reduction_type, reduce_map, reduction_type_, "reduction_type",
              "FrameTransformReductionMgr");
  check_target_downstream_dist_(reduction_type_, target_downstream_dist,
                                component_transform);
  for (int i = 0; i<3; i++){
    component_transform_[i] = component_transform[i];
    target_downstream_dist_[i] = target_downstream_dist[i];
  }

}

//----------------------------------------------------------------------

// directly pupping an enum value results in errors on some systems
template<typename T>
void PUPenum_(PUP::er &p, T &enum_val){
  if (p.isUnpacking()){
    int temp;
    p|temp;
    enum_val =  static_cast<T>(temp);
  } else {
    int temp = static_cast<int>(enum_val);
    p|temp;
  }
}

void FrameTransformReductionMgr::pup(PUP::er &p){
  p|weight_field_;
  p|weight_threshold_;
  PUPenum_<threshold_enum>(p, threshold_type_);
  PUPenum_<frame_trans_reduce_enum>(p, reduction_type_);
  PUParray(p,target_downstream_dist_,3);
  PUParray(p,component_transform_,3);
}

//----------------------------------------------------------------------

void FrameTransformReductionMgr::launch_weighted_average_reduction_
(Block *block, double time_to_next_transform, precision_type precision,
 CkCallback &cb) const throw()
{
  // we refer to the product of weight_field values and cell volumes as
  // "masses" and the product of these "masses" with velocity components as
  // "momentum" components

  // to make this calculation deterministic, going to accumulate masses and
  // momenta of each block before summing them all together

  // Assign the block an index based on it's location in the grid. This
  // will need to be modified once AMR and solvers are in use.
  int ix, iy, iz, nx, ny, nz;
  block->index_global(&ix, &iy, &iz, &nx, &ny, &nz);
  int index = ix + nx*(iy + ny*iz);

  double message_arr[5] = {0., 0., 0., 0., 0.};
  message_arr[0] = (double)index;
  double *mass = &(message_arr[1]);
  double *momentum = &(message_arr[2]);

  auto function =
    [=](double thresh_factor, double weight_val, double v[3],
        int ix, int iy, int iz) {
      // since the cell volume is always constant throughout a block, wait
      // multiply by volume until after processing full block
      double cur_mass = thresh_factor * weight_val;
      (*mass) += cur_mass;
      for (int i = 0; i < 3; i++){
        momentum[i] += v[i] * cur_mass;
      }
    };

  // calc mass and momentum
  if (block->is_leaf()) { local_reduction_(block, function, precision); }

  for (int i=0; i<3; i++) { if (!component_transform_[i]) momentum[i] = 0; }

  // multiply mass and momentum by the volume
  double hx,hy,hz;
  block->data()->field_cell_width(&hx,&hy,&hz);
  double volume = hx*hy*hz;
  (*mass) *= volume;
  for (int j=0; j<3; j++) { momentum[j] *= volume; }

  // now, call the reduction to gather momentum and mass from all blocks
  // (In principle, we could save bandwidth and cpu time by only performing
  // the reduction on the relevant momentum components)
  block->contribute(5*sizeof(double), message_arr, CkReduction::concat,
                    cb);
}

//----------------------------------------------------------------------

void FrameTransformReductionMgr::launch_min_reduction_
(Block *block, double time_to_next_transform, precision_type precision,
 CkCallback &cb) const throw()
{
  // in both cases compute the min on each component
  constexpr double max_dbl = std::numeric_limits<double>::max();
  double velocity[3] = {max_dbl, max_dbl, max_dbl};

  // calculate the minimum velocity components for the local block
  auto function =
    [=, &velocity] (double thresh_factor, double weight_val, double v[3],
                    int ix, int iy, int iz) {
      // weight_val was used to compute thresh_factor)
      double f = thresh_factor;
      for (int i = 0; i < 3; i++){
        velocity[i] = std::fmin(v[i] * f + (1.-f) * max_dbl, velocity[i]);
      }
    };

  if (block->is_leaf()) { local_reduction_(block, function,precision); }

  for (int i=0; i<3; i++) { if (!component_transform_[i]) velocity[i] = 0; }

  // optionally apply a floor of zero
  if (reduction_type_ == frame_trans_reduce_enum::min_zero_floor){
    for (int i = 0; i < 3; i++) {velocity[i] = std::fmax(velocity[i], 0.);}
  }

  // launch the reduction
  block->contribute(3*sizeof(double), velocity, CkReduction::min_double,
                    cb);
}

//----------------------------------------------------------------------


void FrameTransformReductionMgr::launch_target_downstream_dist_reduction_
(Block *block, double time_to_next_transform, precision_type precision,
 CkCallback &cb) const throw()
{

  ASSERT("FrameTransformReductionMgr::min_projected_downstream_dist_",
         "time_to_next_transform must be positive.",
         time_to_next_transform>0);

  // - downstream_dist is the distance from the lower edge of the domain along
  //   a given axis
  // - projected_downstream_dist is the expected downstream distance between
  //   now and the next timestep when the frame velocity will be updated again.
  //   Over this interval, velocities are assumed to be constant.
  const int rank = cello::rank();

  const Config *config_ptr = cello::config();
  double domain_lower[3] = {0., 0., 0.};
  for (int i = 0; i< rank; i++){
    domain_lower[i] = config_ptr->domain_lower[i];
  }

  // Compute widths and offset. These are defined such that the center of cell
  // along axis ax at index i_ax is given by i_ax * widths + offsets;
  double widths[3], offset[3];
  Data* data = block->data();
  data->field_cell_width(&(widths[0]), &(widths[1]), &(widths[2]));
  data->lower(&(offset[0]), &(offset[1]), &(offset[2]));
  for (int i=0; i < 3; i++) {offset[i] += widths[i]*0.5;}

  Field field = block->data()->field();
  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  constexpr double max_dbl = std::numeric_limits<double>::max();
  double proj_downstream_dist[3] = {max_dbl, max_dbl, max_dbl};

  auto function =
    [=,&proj_downstream_dist](double thresh_factor, double weight_val,
                              double v[3], int ix, int iy, int iz)
    {
      // compute the current cell location
      double cur_loc[3];
      cur_loc[0] = ((double)(ix-gx))*widths[0] + offset[0];
      cur_loc[1] = ((double)(iy-gy))*widths[1] + offset[1];
      cur_loc[2] = ((double)(iz-gz))*widths[2] + offset[2];

      double f = thresh_factor;
      for (int i = 0; i < 3; i++){

        double cur_downstream_d = cur_loc[i] - domain_lower[i];
        double cur_proj = cur_downstream_d + v[i] * time_to_next_transform;
        proj_downstream_dist[i] =
          std::fmin(proj_downstream_dist[i],
                    cur_proj * f + (1.-f) * max_dbl);
      }
    };

  double velocity[3] = {max_dbl, max_dbl, max_dbl};
  if (block->is_leaf()) {
    local_reduction_(block, function, precision);
    for (int i =0; i < 3; i++){
      if (proj_downstream_dist[i] < max_dbl){
        velocity[i] = ((proj_downstream_dist[i] - target_downstream_dist_[i])
                       / time_to_next_transform);
      }
    }
  }

  block->contribute(3*sizeof(double), velocity, CkReduction::min_double, cb);

}

//----------------------------------------------------------------------

void FrameTransformReductionMgr::extract_final_velocity
(CkReductionMsg * msg, double (&v)[3], bool (&update_component)[3])
  const throw()
{

  for (int i = 0; i<3; i++){
    v[i] = 0.;
    update_component[i] = ( (i < cello::rank()) &&
                            component_transform_[i] );
  }
  double *values=(double *)msg->getData();

  if (reduction_type_ == frame_trans_reduce_enum::weighted_average){

    // To make our results deterministic, extract all the concatenated
    // reduction messages and sort them by the block index. We sort in place
    // to avoid memory allocation (the charm++ ampi implementation does this
    // too)

    int n = msg->getSize()/(5*sizeof(double));

    for (int i = 0; i < n; i++){
      int index = (int) values[i*5]; // this had been cast to double for
                                     // convenience
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
    for (int i = 0; i<3; i++){
      if (mass > 0){
        v[i] = momentum[i]/mass;
      } else {
        update_component[i] = false;
      }
    }

  } else {
    ASSERT("FrameTransformReductionMgr::extract_final_velocity",
           "CkReductionMsg is expected to hold 3 doubles",
           msg->getSize() == 3 *sizeof(double));

    for (int i = 0; i < 3; i++) {
      // if no cells met the threshold, the values will be equal to
      // std::numeric_limits<double>::max()
      double max_dbl = std::numeric_limits<double>::max();
      if (values[i] != max_dbl){
        v[i] = values[i];
      } else {
        update_component[i] = false;
      }
    }
  }
}

//----------------------------------------------------------------------
template <class Function>
void FrameTransformReductionMgr::local_reduction_
(Block * block, Function func, precision_type precision) const throw()
{
  switch (precision) {
  case precision_single:
    { local_reduction_helper_<float,Function>(block, func); }
    break;
  case precision_double:
    { local_reduction_helper_<double,Function>(block, func); }
    break;
  case precision_extended80:
  case precision_extended96:
  case precision_quadruple:
    { local_reduction_helper_<long double,Function>(block, func); }
    break;
  }
}

template <class T, class Function>
void FrameTransformReductionMgr::local_reduction_helper_
(Block * block, Function func) const throw()
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

  T* weight_vals = (T *) field.values(weight_field_);
  T * v3[3] = { (T*) (              field.values("velocity_x")),
                (T*) ((rank >= 2) ? field.values("velocity_y") : NULL),
                (T*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  // define lambda function to apply threshold.
  const double weight_threshold = weight_threshold_;
  const threshold_enum threshold_type = threshold_type_;

  auto satisfies_thresh =
    [=](double value)->bool{
      switch(threshold_type) {
      case(threshold_enum::ignore)      : {return true;}
      case(threshold_enum::lower_limit) : {return value >= weight_threshold;}
      case(threshold_enum::upper_limit) : {return value <= weight_threshold;}
      }
      return false;
    };

  // only contain material in active zone
  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {

        int i = ix + ndx*(iy + ndy*(iz));

        double weight_val = (double) weight_vals[i];
        double thresh_factor = satisfies_thresh(weight_val);
        double cur_v[3] = {0., 0., 0.};
        for (int j=0; j<rank; j++){ cur_v[j] = (double) v3[j][i]; }

        func(thresh_factor, weight_val, cur_v, ix, iy, iz);
      }
    }
  }
}

//======================================================================

MethodFrameTransform::MethodFrameTransform
(const bool (&component_transform)[3], std::string weight_field,
 double weight_threshold, std::string threshold_type,
 std::string reduction_type, const double (&target_downstream_dist)[3])
  : Method(),
    reduction_mgr_(weight_field, weight_threshold, threshold_type,
                   reduction_type, target_downstream_dist, component_transform)
{

  // Copy component_transform entries
  int num_components = 0;
  for (int i=0; i<3; i++){
    num_components += (int)(component_transform[i]);
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
	   !(component_transform[1] || component_transform[2]));
  } else if (rank == 2){
    ASSERT("MethodFrameTransform",
	   ("The z velocity components can't be transformed for a 2D"
	    "simulation."), !component_transform[2]);
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
}

//----------------------------------------------------------------------

void MethodFrameTransform::pup(PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  p|reduction_mgr_;
}

//----------------------------------------------------------------------

void MethodFrameTransform::compute( Block * block) throw()
{
  TRACE_FRAME_TRANSFORM;

  ASSERT("MethodFrameTransform::compute", "Not currently compatible with AMR",
	 block->level() == 0 && block->is_leaf());

  // compute the time to the next invocation of MethodFrameTransform::compute
  double time_to_next_transform = -1;
   if ((passive_schedule_ != nullptr) &&
       ( (passive_schedule_->type() == schedule_type_time) ||
         (passive_schedule_->type() == schedule_type_minimum_time) ) ){
     double time_at_cycle_end = block->time() + block->dt();
     time_to_next_transform =
       passive_schedule_->time_next() - time_at_cycle_end;
  }

  Field field = block->data()->field();

  // construct the charm++ callback to indicate that
  // p_method_frame_transform_end should be called after the reduction
  CkCallback cb(CkIndex_Block::p_method_frame_transform_end(NULL),
		block->proxy_array());

  reduction_mgr_.launch_reduction(block, time_to_next_transform,
                                  field_precision_(field), cb);
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

  // determine the new frame velocity (measured in the current frame) from the
  // the result of the reduction message
  double frame_velocity[3] = {0., 0., 0.};
  bool update_component[3] = {false, false, false};
  reduction_mgr_.extract_final_velocity(msg, frame_velocity,
                                        update_component);

  // Update the velocity of the current frame measured with respect to the
  // frame when the gas was originally initialized (add the new frame velocity
  // to the old frame velocity)
  double vx,vy,vz;
  block->frame_velocity(&vx, &vy, &vz);

  double new_v[3];
  new_v[0] = (update_component[0]) ? vx + frame_velocity[0] : vx;
  new_v[1] = (update_component[1]) ? vy + frame_velocity[1] : vy;
  new_v[2] = (update_component[2]) ? vz + frame_velocity[2] : vz;

  // the update time is recorded as the time+dt becasue the update occurs at
  // the end of the cycle (all other preceeding methods were applied using the
  // preceeding velocity). Another benefit of this choice is that if an output
  // is written at the end of this cycle, both the time of the last velocity
  // update AND the last updated origin offset will reflect the actual
  // simulation time AND the actual origin offset.
  block->update_frame_properties(block->time() + block->dt(),
                                 new_v[0], new_v[1], new_v[2]);

  // summarize frame properties (only print if using the root block)
  if (block->index().is_root()) {
    double offset[3];
    block->last_updated_origin_offset(&(offset[0]), &(offset[1]),
                                      &(offset[2]));

    const char* axes = "xyz";
    const char* prefixes[] = {"   Post", "Skipped"};
    for (int i = 0; i < cello::rank(); i++){

      const char* prefix = (update_component[i]) ? prefixes[0] : prefixes[1];
      cello::monitor()->print
        ("Method",
         "Ref-Frame %s-%c-transform: vel = %+20.16e origin-offset = %20.16e",
         prefix, axes[i], new_v[i], offset[i]);
    }
  }

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

    // We are iterating over the entire mesh. To just iterate over active zone,
    // don't modify nx,ny,nz and set ix0, iy0, iz0 to gx, gy, gz
    int ix0 = 0, iy0 = 0, iz0 = 0;

    nx = ndx;
    ny = ndy;
    nz = ndz;

    for (int i = 0; i < 3; i++){
      if (!update_component[i]) {
        frame_velocity[i] = 0;
      }
    }

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

    if (( (i <  3) && reduction_mgr_.component_to_be_transformed(i)  ) ||
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


double MethodFrameTransform::timestep ( Block * block ) const throw()
{
  double out = std::numeric_limits<double>::max();

  // this is a somewhat hacky workaround. We would be better off updating
  // the control_stopping code
  if ( (passive_schedule_ != nullptr) &&
       (passive_schedule_->type() == schedule_type_time) ){
    ASSERT("MethodFrameTransform::timestep",
           "courant must be 1", courant_ == 1.);
    ASSERT("MethodFrameTransform::timestep",
           "courant_global must be 1",
           Method::courant_global == 1.);

    // we need to advance the schedule to the current time and cycle
    // it's ok to do this multiple times in a row as long as the arguments don't
    // change between calls
    passive_schedule_->is_scheduled(block->cycle(), block->time());

    out = passive_schedule_->time_next() - block->time();
    // SANITY CHECK
    ASSERT("MethodFrameTransform::timestep", "The timestep can't be negative.",
           out >= 0.);
  }

  return out;
}

