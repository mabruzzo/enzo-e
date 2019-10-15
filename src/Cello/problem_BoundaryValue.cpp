// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-02
/// @brief    Implementation of the default BoundaryValue boundary value class

#include "problem.hpp"

//----------------------------------------------------------------------

void update_boundary_mask_(bool mask[3][2], axis_enum axis, face_enum face,
			   bool value){
  if (axis==axis_all && face==face_all) {
    for (int axis_ind=0; axis_ind<3; axis_ind++) {
      for (int face_ind=0; face_ind<2; face_ind++) {
        mask[axis_ind][face_ind] = value;
      }
    }
  } else if (axis==axis_all) {
    for (int axis_ind=0; axis_ind<3; axis_ind++) {
      mask[axis_ind][face] = value;
    }
  } else if (face==face_all) {
    for (int face_ind=0; face_ind<2; face_ind++) {
      mask[axis][face_ind] = value;
    }
  } else {
    mask[axis][face] = value;
  }
}

//----------------------------------------------------------------------

BoundaryValue::BoundaryValue(axis_enum axis, face_enum face, Value * value, 
			     std::vector<std::string> field_list,
			     bool possible_velocity_frame_transform,
			     const std::vector<Boundary*> * earlier_boundaries)
  throw()
  : Boundary(axis,face,0),
    value_(value),
    field_list_(field_list)
{
  // Initialize velocity_frame_transform_ values
  for (int axis_ind=0; axis_ind<3; axis_ind++) {
    for (int face_ind=0; face_ind<2; face_ind++) {
      velocity_frame_transform_[axis_ind][face_ind] = false;
    }
  }

  // If not there is no necessary velocity frame transformation, return now
  if (!possible_velocity_frame_transform) { return; }

  // Check that this instance enforces a boundary for one of the fields that
  // must be transformed ("velocity_x", "velocity_y", "velocity_z" or
  // "total_energy")
  std::string check_fields[4] = {"velocity_x", "velocity_y", "velocity_z",
				 "total_energy"};
  bool has_field = false;
  for (int i=0; i < 4; i++){
    if (field_list_.end() != std::find(field_list_.begin(), field_list_.end(),
				       check_fields[i])){ 
      has_field = true;
    }
  }

  // If this instance doesn't handle one of the above fields return now
  if (!has_field) {return;}

  // Update the appropriate values of velocity_frame_transform_ to indicate
  // the boundary faces that this instance handles transformations for
  update_boundary_mask_(velocity_frame_transform_, axis, face, true);

  // Finally go through any Boundary instances that get enforced before this
  // instance and make sure no other BoundaryValue instances applies the
  // velocity transform for the same boundary faces as this instance.
  ensure_exclusive_transform_(earlier_boundaries);
}

//----------------------------------------------------------------------

// Iterates over earlier_boundaries (a vector of pointers to Boundary
// objects that are enforced before *this) and updates the attributes of
// any BoundaryValue instance to ensure that they do not apply the velocity
// transform for any of the boundaries that *this is responsible for.
void BoundaryValue::ensure_exclusive_transform_(const std::vector<Boundary*> *
						earlier_boundaries)
  throw ()
{
  std::string ref_name = name();

  for (std::size_t i = 0; i<earlier_boundaries->size(); i++){
    Boundary * boundary_ptr = (*earlier_boundaries)[i];

    // Skip the Boundary condition if it is not an instance of BoundaryValue
    // or it does not apply to any of the boundary faces that *this applies to
    if (boundary_ptr->name() != ref_name) {continue; }
    BoundaryValue* comp_boundary = dynamic_cast<BoundaryValue*>(boundary_ptr);
    face_enum ptr_face = comp_boundary->face_;
    if ((ptr_face != face_all) && (face_ != face_all) && (ptr_face != face_)){
      continue;
    }
    axis_enum ptr_axis = comp_boundary->axis_;
    if ((ptr_axis != axis_all) && (axis_ != axis_all) && (ptr_axis != axis_)){
      continue;
    }

    update_boundary_mask_(comp_boundary->velocity_frame_transform_, axis_,
			  face_, false);    
  }
}

//----------------------------------------------------------------------

void BoundaryValue::enforce 
(Block * block, face_enum face, axis_enum axis) const throw()
{
  if ( ! applies_(axis,face)) return;

  if (face == face_all) {
    enforce(block,face_lower,axis);
    enforce(block,face_upper,axis);
  } else if (axis == axis_all) {
    enforce(block,face,axis_x);
    enforce(block,face,axis_y);
    enforce(block,face,axis_z);
  } else {

    Data * data = block->data();
    Field field = data->field();

    if ( ! field.ghosts_allocated() ) {
      ERROR("EnzoBoundary::enforce",
	    "Function called with ghosts not allocated");
    }

    double xm,ym,zm;
    double xp,yp,zp;
    data -> lower(&xm,&ym,&zm);
    data -> upper(&xp,&yp,&zp);

    double t = block->time();

    for (size_t index = 0; index < field_list_.size(); index++) {

      int nx,ny,nz;
      field.size(&nx,&ny,&nz);

      int index_field = field.field_id(field_list_[index]);
      int gx,gy,gz;
      field.ghost_depth(index_field,&gx,&gy,&gz);

      int cx,cy,cz;
      field.centering(index_field, &cx,&cy,&cz);

      int ndx=nx+2*gx+cx;
      int ndy=ny+2*gy+cy;
      int ndz=nz+2*gz+cz;

      double * x = new double [ndx];
      double * y = new double [ndy];
      double * z = new double [ndz];

      data->field_cell_faces(x,y,z,gx,gy,gz,cx,cy,cz);

      void * array = field.values(index_field);

      precision_type precision = field.precision(index_field);

      int ix0=0, iy0=0, iz0=0;

      nx = ndx;
      ny = ndy;
      nz = ndz;

      if (axis == axis_x) nx=gx;
      if (axis == axis_y) ny=gy;
      if (axis == axis_z) nz=gz;

      if (face == face_upper) {
	if (axis == axis_x) ix0 = ndx - gx;
	if (axis == axis_y) iy0 = ndy - gy;
	if (axis == axis_z) iz0 = ndz - gz;
      }

      int i0=ix0 + ndx*(iy0 + ndy*iz0);

      bool * mask = 0;

      if (mask_ != nullptr) mask = new bool [nx*ny*nz];

      switch (precision) {
      case precision_single:
	{
	  float * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (float *)array;
	    array = new float [ndx*ndy*ndz];
	  }
	  
	  value_->evaluate((float *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((float *)temp)[i]=((float *)array)[i];
	    delete [] ((float*)array);
	    array = temp;
	  }
	}
       	break;
      case precision_double:
	{
	  double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (double *)array;
	    array = new double [ndx*ndy*ndz];
	  }
	  value_->evaluate((double *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) ((double *)temp)[i]=((double *)array)[i];
	    delete [] ((double *)array);
	    array = temp;
	  }
	}
       	break;
      case precision_extended80:
      case precision_extended96:
      case precision_quadruple:
	{
	  long double * temp = 0;
	  if (mask_ != nullptr) {
	    temp = (long double *)array;
	    array = new long double [ndx*ndy*ndz];
	  }
	  value_->evaluate((long double *)array+i0, t, 
			   ndx,nx,x+ix0, 
			   ndy,ny,y+iy0,
			   ndz,nz,z+iz0);
	  if (mask_ != nullptr) {
	    for (int i=0; i<ndx*ndy*ndz; i++) 
	      ((long double *)temp)[i]=((long double *)array)[i];
	    delete [] ((long double *)array);
	    array = temp;
	  }
	}
       	break;
      }

      delete [] x;
      delete [] y;
      delete [] z;
      delete [] mask;
    }
    if (velocity_frame_transform_[axis][face]) {
      // update the reference frame of the velocity
      transform_velocity_energy_(block, face, axis);
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void BoundaryValue::copy_(T * field, double * value,
			  int ndx, int ndy, int ndz,
			  int nx,  int ny,  int nz,
			  int ix0, int iy0, int iz0) const throw()
{
  for (int ix=ix0; ix<ix0+nx; ix++) {
    for (int iy=iy0; iy<iy0+ny; iy++) {
      for (int iz=iz0; iz<iz0+nz; iz++) {
	int iv = (ix-ix0) + nx*((iy-iy0) + ny*(iz-iz0));
	int ib = ix + ndx*(iy + ndy*(iz));
	field[ib] = (T) value[iv];
      }
    }
  }
}

//----------------------------------------------------------------------


void BoundaryValue::transform_velocity_energy_(Block * block, face_enum face,
					       axis_enum axis)
  const throw()
{
  // this method can only be called if velocity_x is a valid field
  Field field = block->data()->field();

  double frame_velocity[3];
  block->frame_velocity(frame_velocity, frame_velocity+1, frame_velocity+2);

  if ((frame_velocity[0] == 0.) && (frame_velocity[1] == 0.) &&
      (frame_velocity[2] == 0.)){
    return;
  }

  std::string field_names[4] = {"velocity_x", "velocity_y", "velocity_z",
				"total_energy"};

  precision_type precision = field.precision(field.field_id("velocity_x"));

  
  for (int i=0; i<4; i++){
    std::string field_name = field_names[i];
    if (!field.is_field(field_name)) { continue; }

    int index_field = field.field_id(field_name);
    ASSERT1("BoundaryValue",
	    ("\"velocity_x\" and \"%s\" must have the same precision for "
	     "reference frame updates."),
	    precision == field.precision(index_field), field_name.c_str());
    ASSERT1("BoundaryValue",
	    "\"%s\" must be cell-centered for reference frame updates.",
	    field.is_centered(index_field), field_name.c_str());
  }

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int ndx=nx+2*gx;
  int ndy=ny+2*gy;
  int ndz=nz+2*gz;

  int ix0=0, iy0=0, iz0=0;

  nx = ndx;
  ny = ndy;
  nz = ndz;

  if (axis == axis_x) nx=gx;
  if (axis == axis_y) ny=gy;
  if (axis == axis_z) nz=gz;

  if (face == face_upper) {
    if (axis == axis_x) ix0 = ndx - gx;
    if (axis == axis_y) iy0 = ndy - gy;
    if (axis == axis_z) iz0 = ndz - gz;
  }

  switch (precision) {
  case precision_single:
    {
      float vx = (float) frame_velocity[0];
      float vy = (float) frame_velocity[1];
      float vz = (float) frame_velocity[2];
      MethodFrameTransform::transform_field (block, vx, vy, vz,
					     ndx, ndy, ndz,
					     nx,  ny,  nz,
					     ix0, iy0, iz0,
					     0);
      break;
    }
  case precision_double:
    {
      double vx = frame_velocity[0];
      double vy = frame_velocity[1];
      double vz = frame_velocity[2];
      MethodFrameTransform::transform_field (block, vx, vy, vz,
					     ndx, ndy, ndz,
					     nx,  ny,  nz,
					     ix0, iy0, iz0,
					     0);
      break;
    }
  case precision_extended80:
  case precision_extended96:
  case precision_quadruple:
    {
      long double vx = (long double) frame_velocity[0];
      long double vy = (long double) frame_velocity[1];
      long double vz = (long double) frame_velocity[2];
      MethodFrameTransform::transform_field (block, vx, vy, vz,
					     ndx, ndy, ndz,
					     nx,  ny,  nz,
					     ix0, iy0, iz0,
					     0);
      break;
    }
  }
}
