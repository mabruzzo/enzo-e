// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodFrameTransform.hpp
/// @author   Matthew W. Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-09-23
/// @brief    [\ref Enzo]  Declares the MethodFrameTransform class

#ifndef PROBLEM_METHOD_FRAME_TRANSFORM
#define PROBLEM_METHOD_FRAME_TRANSFORM

/// @enum     threshold_enum
/// @brief    Describes how to treat weight_threshold
enum threshold_enum {
  ignore = 0,
  lower_limit,
  upper_limit
};

/// @enum     frame_trans_reduce_enum
/// @brief    Describes the type of reduction
enum frame_trans_reduce_enum {
  unknown = 0,
  weighted_average,
  min,
  min_zero_floor,
};

class FrameTransformReductionMgr{
  // This is a helper class that manages operations associated with the block
  // reduction for calculating incremental changes in the frame velocity. In
  // all cases, weight_field is used with the threshold information to identify
  // which cells to involve in the calculation 
  //
  // The currently supported reduction types are:
  //     - "weighted_average": each component of the new frame velocity is the
  //       average of the velocity component of all selected cells, weighted
  //       by weight_field. Note: weight_field is assumed to be a density-like
  //       quantity (it will be multplied by cell-volume)
  //     - "min": each component of the new frame velocity is the minimum value
  //       of the velocity component of all selected cells
  //     - "min_zero_floor": same as "min", except the velocity is not allowed
  //       a floor of zero is placed on the final value

public:
  /// Main Constructor
  FrameTransformReductionMgr(std::string weight_field, double weight_threshold,
                             std::string threshold_type,
                             std::string reduction_type) throw();

  /// Default Constructor - mainly used for charm++ unpacking
  FrameTransformReductionMgr() throw()
    : weight_field_(""), weight_threshold_(0.0),
      threshold_type_(threshold_enum::ignore),
      reduction_type_(frame_trans_reduce_enum::unknown)
  { }

  /// CHARM++ Pack / Unpack function
  ///
  /// Since the class is not meant to be referred to by pointer, it doesn't
  /// need to be declared PUPable
  void pup(PUP::er &p);

  /// Launch the charm++ reduction to compute the new frame velocity
  ///
  /// @tparam T the precision of the block fields
  template<typename T>
  void launch_reduction(Block * block, const bool component_transform[3],
                        CkCallback &cb) const throw();

  /// Function called after completion of the reduction to determine the new
  /// frame velocity (in the current reference frame) from the reduction result
  void extract_final_velocity(CkReductionMsg * msg, double v[3]) const throw();

private: // helper functions

  /// helper function used to compute the reduction values for the local block
  template <class T, class Function>
  void local_block_reduction_(Block * block, Function func) const throw();

private: // attributes
  /// Name of the density field used for weighting the velocity.
  std::string weight_field_;

  /// A threshold for the weight_field
  double weight_threshold_;

  /// How the threshold should be applied
  threshold_enum threshold_type_;

  /// The type of reduction to be applied
  frame_trans_reduce_enum reduction_type_;
};
  

class MethodFrameTransform : public Method {

  /// @class    MethodFrameTransform
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Transforms the reference frame such that the
  ///           designated components of the weighted average velocity of
  ///           are zero. The average is weighted by some kind of density field
  ///           (whether its density, a passive scalar, or some other field
  ///           measuring a quantity per unit volume)

public: // interface
  /// Create a new MethodFrameTransform object
  MethodFrameTransform(bool component_transform[3],
		       std::string weight_field,
		       double weight_threshold,
		       std::string threshold_type,
                       std::string reduction_type);

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodFrameTransform);

  /// Charm++ PUP::able migration constructor
  MethodFrameTransform (CkMigrateMessage *m)
    : Method (m),
      component_transform_(), // initialize to (false,false,false)
      reduction_mgr_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "frame_transform"; }

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

protected: // attributes

  /// Indicates the velocity components to apply the transformation for
  /// - the corresponding components are ordered (x,y,z)
  bool component_transform_[3];

  /// Manages the reduction to compute incremental changed in frame velocity
  FrameTransformReductionMgr reduction_mgr_;
};

#endif /* PROBLEM_METHOD_FRAME_TRANSFORM */
