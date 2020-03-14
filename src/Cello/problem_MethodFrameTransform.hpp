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
		       bool cycle_based_update,
		       double update_start, double update_step,
		       double weight_threshold,
		       std::string threshold_type);

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodFrameTransform);

  /// Charm++ PUP::able migration constructor
  MethodFrameTransform (CkMigrateMessage *m)
    : Method (m),
      component_transform_(), // initialize to (false,false,false)
      weight_field_(""),
      weight_threshold_(0.0),
      threshold_type_(threshold_enum::ignore),
      cycle_based_update_(true),
      cycle_start_(0),
      cycle_step_(1),
      time_start_(0.0),
      time_step_(0.0),
      prev_time_update_(0.)
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

  /// Checks whether current simulation time/cycle preceedes the first
  /// scheduled frame velocity update
  bool precede_first_transform_(Block * block) const throw();

  /// Checks if a frame velocity update is scheduled for the current simulation
  /// time/cycle preceedes. If an update is scheduled and scheduling is
  /// time-based, this records that the current simulation time as the time at
  /// which the most recent velocity update was performed.
  ///
  /// @details
  /// Because MethodFrameTransform can't rely on the method being executed at
  /// specific simulation times, time scheduled frame velocity updates are
  /// performed at the simulation times immediately following the scheduled
  /// times. Consequently, the simulation time at which the most recent frame
  /// velocity update was performed must be recorded (to avoid multiple updates
  /// for a single scheduled time). This function effectively assumes that
  /// the frame velocity will be updated for any cycle that it returns true,
  /// which is why it records that the current simulation time as the time at
  /// which the most recent velocity update was performed.
  bool scheduled_velocity_update_(Block * block) throw();

  /// Returns the precision of the velocity and energy fields. Makes sure that
  /// the precision is the same for each of them.
  precision_type field_precision_(Field &field) const throw();

  /// Computes the sum of the weight_field multiplied by cell volume ("mass")
  /// and the "mass" multiplied by the velocity ("momementum") over the block.
  template <class T>
  void block_totals_(Block *block, double &mass,
		     double momentum[3]) const throw();

protected: // attributes

  /// Indicates the velocity components to apply the transformation for
  /// - the corresponding components are ordered (x,y,z)
  bool component_transform_[3];

  /// Name of the density field used for weighting the velocity.
  std::string weight_field_;

  /// A threshold for the weight_field
  double weight_threshold_;

  /// How the threshold should be applied
  threshold_enum threshold_type_;

  /// Whether updates are based on cycle OR simulation time
  bool cycle_based_update_;

  /// Cycle to start tracking the velocity and perform the transform
  int cycle_start_;
  /// the number of cycles to wait between updating the reference frame
  int cycle_step_;

  /// simulation time to start tracking the velocity and perform the transform
  double time_start_;
  /// the amount of simulation time between updating the reference frame
  double time_step_;
  /// Records the time at which the previous velocity update was made
  double prev_time_update_;

};

#endif /* PROBLEM_METHOD_FRAME_TRANSFORM */
