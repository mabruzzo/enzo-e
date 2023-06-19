// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlong.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the EnzoProlong class
///
/// This class serves to encapsulate Enzo's interpolate() function

#ifndef PROBLEM_ENZO_PROLONG_HPP
#define PROBLEM_ENZO_PROLONG_HPP

class EnzoProlong : public Prolong {

  /// @class    EnzoProlong
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  EnzoProlong(std::string method,
              bool positive,
              bool use_linear) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlong);

  /// CHARM++ migration constructor
  EnzoProlong(CkMigrateMessage *m)
    : Prolong(m),method_(-1),
      positive_(1),
      use_linear_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Prolong fine Field values in the child block (icx,icy,icz) to parent
  virtual void apply 
  ( precision_type precision,
    void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false);

  /// Return the name identifying the prolongation operator
  virtual std::string name () const { return "enzo"; }

protected: // protected virtual methods
  
  /// Amount of padding required in coarse region (default 0)
  virtual int coarse_padding_() const
  { return 1; }

private: // functions

  void apply_
  ( enzo_float *  values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const enzo_float * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false);
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Interpolation Method: see Enzo documenation
  int method_;

  /// Positivity flag
  int positive_;

  /// Whether to bypass and actually use ProlongLinear (for debugging)
  bool use_linear_;

};

#endif /* PROBLEM_ENZO_PROLONG_HPP */

