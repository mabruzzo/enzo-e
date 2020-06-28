// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-04-01
/// @brief    [\ref Problem] Declaration for the BoundaryValue component

#ifndef PROBLEM_BOUNDARY_VALUE_HPP
#define PROBLEM_BOUNDARY_VALUE_HPP

class BoundaryValue : public Boundary
{

  /// @class    BoundaryValue
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Encapsulate a BoundaryValue conditions generator

public: // interface

  /// Create a new BoundaryValue
  BoundaryValue() throw() 
  : Boundary (), value_(NULL), field_list_()
  {
    for (int axis_ind=0; axis_ind<3; axis_ind++) {
      for (int face_ind=0; face_ind<2; face_ind++) {
	velocity_frame_transform_[axis_ind][face_ind] = false;
      }
    }
  }

  /// Create a new BoundaryValue
  BoundaryValue(axis_enum axis, face_enum face, Value * value, 
		std::vector<std::string> field_list,
		bool possible_velocity_frame_transform = false,
		const std::vector<Boundary*> * earlier_boundaries = 0) throw();

  /// Destructor
  virtual ~BoundaryValue() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(BoundaryValue);

  BoundaryValue(CkMigrateMessage *m)
    : Boundary (m),
      value_(NULL),
      field_list_()
  {
    for (int axis_ind=0; axis_ind<3; axis_ind++) {
      for (int face_ind=0; face_ind<2; face_ind++) {
	velocity_frame_transform_[axis_ind][face_ind] = false;
      }
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    // NOTE: change this function whenever attributes change
    Boundary::pup(p); 
    TRACEPUP; 

    p | *value_;
    p | field_list_;

    for (int axis_ind=0; axis_ind<3; axis_ind++){
      PUParray(p,velocity_frame_transform_[axis_ind],2);
    }
 
  };

public: // virtual functions

  /// Enforce BoundaryValue conditions

  virtual void enforce (Block   * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw();

  /// Return the name of this boundary
  virtual std::string name () throw()
  { return "value"; }

protected: // functions

  template <class T>
  void copy_(T * field, double * value,
	     int ndx, int ndy, int ndz,
	     int nx,  int ny,  int nz,
	     int ix0, int iy0, int iz0) const throw ();

  /// Ensures that `this` is the only BoundaryValue instance to apply the
  /// velocity reference frame transformation along its specified boundary
  /// face(s) by modifying relevant (previously initialized) instances of
  /// BoundaryValue. This is called during initialization
  void ensure_exclusive_transform_(const std::vector<Boundary*> *
				   earlier_boundaries) throw ();

  void transform_velocity_energy_(Block * block, face_enum face,
				  axis_enum axis) const throw();

protected: // attributes

  Value * value_;
  std::vector<std::string> field_list_;

  /// If the velocity frame tranform should be applied on a given face
  bool velocity_frame_transform_[3][2];

};

#endif /* PROBLEM_BOUNDARY_VALUE_HPP */
