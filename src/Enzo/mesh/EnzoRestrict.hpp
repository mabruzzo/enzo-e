// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRestrict.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    [\ref Problem] Declaration of the EnzoRestrict class
///
/// This class serves to encapsulate Enzo's restriction operations

#ifndef PROBLEM_ENZO_RESTRICT_HPP
#define PROBLEM_ENZO_RESTRICT_HPP

class EnzoRestrict : public Restrict {
  /// @class    EnzoRestrict
  /// @ingroup  Problem
  /// @brief    [\ref Problem]

public:  // interface
  /// Constructor
  EnzoRestrict(std::string method) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoRestrict);

  /// CHARM++ migration constructor
  EnzoRestrict(CkMigrateMessage* m) : Restrict(m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p);

  /// Restrict coarse Field values to the child block (icx,icy,icz)

  int apply(precision_type precision, void* values_c, int nd3_c[3],
            int im3_c[3], int n3_c[3], const void* values_f, int nd3_f[3],
            int im3_f[3], int n3_f[3], bool accumulate = false);

  /// Return the name identifying the restrict operator
  virtual std::string name() const { return "enzo"; }

private:  // functions
  int apply_(enzo_float* values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
             const enzo_float* values_f, int nd3_f[3], int im3_f[3],
             int n3_f[3], bool accumulate = false);

private:  // attributes
  // NOTE: change pup() function whenever attributes change
};

#endif /* PROBLEM_ENZO_RESTRICT_HPP */
