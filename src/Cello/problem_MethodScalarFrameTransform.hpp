// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodScalarFrameTransform.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-09-23
/// @brief    [\ref Enzo]  Declares the MethodScalarFrameTransform class

#ifndef PROBLEM_METHOD_SCALAR_FRAME_TRANSFORM
#define PROBLEM_METHOD_SCALAR_FRAME_TRANSFORM

class MethodScalarFrameTransform : public Method {

  /// @class    MethodScalarFrameTransform
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Transforms the reference frame such that the
  ///           designated components of the mass-weighted average velocity of
  ///           a designated passive scalar are zero.

public: // interface
  /// Create a new MethodScalarFrameTransform object
  MethodScalarFrameTransform(bool component_transform[3],
			     std::string passive_scalar,
			     int initial_cycle, int update_stride);

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodScalarFrameTransform);

  /// Charm++ PUP::able migration constructor
  MethodScalarFrameTransform (CkMigrateMessage *m)
    : Method (m),
      component_transform_(), // initialize to (false,false,false)
      passive_scalar_(""),
      initial_cycle_(0),
      update_stride_(1)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "scalar_frame_transform"; }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw();

  /// Performs the reference frame transformation on the fields given the
  /// velocity of the new frame measured with respect to the frame in which the
  /// fields are initially provided
  template <class T>
  static void transform_field(Block * block,
			      T vx,    T vy,    T vz,
			      int ndx, int ndy, int ndz,
			      int nx,  int ny,  int nz,
			      int ix0, int iy0, int iz0,
			      int i_hist) throw();

private: // methods

  /// Returns the precision of the velocity and energy fields. Makes sure that
  /// the precision is the same for each of them.
  precision_type field_precision_(Field &field) const throw();

  /// Sums the total momentum and mass of a passive scalar over the block
  template <class T>
  void block_totals_(Block *block, double &scalar_mass,
		     double scalar_momentum[3]) const throw();

protected: // attributes

  /// Indicates the velocity components to apply the transformation for
  /// - the corresponding components are ordered (x,y,z)
  bool component_transform_[3];

  /// Name of the passively advected scalar to transform
  std::string passive_scalar_;

  /// Cycle to start tracking the velocity and perform the transform
  int initial_cycle_;
  /// the number of cycles to wait between updating the reference frame
  int update_stride_;

};

#endif /* PROBLEM_METHOD_SCALAR_FRAME_TRANSFORM */
