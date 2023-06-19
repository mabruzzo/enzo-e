// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCicInterp.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2016-05-05
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeCicInterp functions

#ifndef ENZO_ENZO_COMPUTE_CIC_INTERP_HPP
#define ENZO_ENZO_COMPUTE_CIC_INTERP_HPP

class EnzoComputeCicInterp : public Compute {

  /// @class    EnzoComputeCicInterp
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate CIC (Cloud-in-cell) particle-field interpolation

public: // interface

  /// Create a new EnzoComputeCicInterp object
  EnzoComputeCicInterp (std::string field_name,
			std::string particle_type,
			std::string particle_attribute,
			double dt = 0.0);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeCicInterp);
  
  /// Charm++ PUP::able migration constructor
  EnzoComputeCicInterp (CkMigrateMessage *m)
    : Compute(m),
      it_p_(0),
      ia_p_(0),
      if_(0),
      dt_(0.0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Perform the computation on the block
  virtual void compute( Block * block) throw();

private: // functions

  void compute_(Block * block);

private: // attributes

  /// particle type
  int it_p_;

  /// particle attribute
  int ia_p_;

  /// id of field type
  int if_;

  /// dt at which to apply the interpolation
  double dt_;

};

#endif /* ENZO_ENZO_COMPUTE_CIC_INTERP_HPP */
