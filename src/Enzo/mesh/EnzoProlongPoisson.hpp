// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlongPoisson.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-01-16
/// @brief    [\ref Enzo] Declaration of the EnzoProlongPoisson class

#ifndef ENZO_PROLONG_POISSON_HPP
#define ENZO_PROLONG_POISSON_HPP

class EnzoProlongPoisson : public Prolong

{
  /// @class    EnzoProlongPoisson
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo]

public:  // interface
  /// Constructor
  EnzoProlongPoisson() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlongPoisson);

  /// CHARM++ migration constructor
  EnzoProlongPoisson(CkMigrateMessage* m) : Prolong(m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p) {
    TRACEPUP;
    Prolong::pup(p);
  }

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent

  virtual void apply(precision_type precision, void* values_f, int nd3_f[3],
                     int im3_f[3], int n3_f[3], const void* values_c,
                     int nd3_c[3], int im3_c[3], int n3_c[3],
                     bool accumulate = false);

private:  // functions
  void apply_(enzo_float* values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
              const enzo_float* values_c, int nd3_c[3], int im3_c[3],
              int n3_c[3], bool accumulate = false);

  /// Return the name identifying the prolongation operator
  virtual std::string name() const { return "poisson"; }

private:  // attributes
  // NOTE: change pup() function whenever attributes change
};

#endif /* ENZO_PROLONG_POISSON_HPP */
