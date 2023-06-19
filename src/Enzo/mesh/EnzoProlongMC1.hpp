// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlongMC1.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    [\ref Problem] Declaration of the EnzoProlongMC1 class
///
/// This class serves to encapsulate Enzo's interpolate() function

#ifndef PROBLEM_ENZO_PROLONG_MC1_HPP
#define PROBLEM_ENZO_PROLONG_MC1_HPP

class EnzoProlongMC1 : public Prolong {

  /// @class    EnzoProlongMC1
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  EnzoProlongMC1(std::string method) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProlongMC1);

  /// CHARM++ migration constructor
  EnzoProlongMC1(CkMigrateMessage *m)
    : Prolong(m), method_(0)
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
  virtual std::string name () const { return "mc1"; }

private: // functions

  void apply_
  ( enzo_float *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
    const enzo_float * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
    bool accumulate = false);
  
private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Interpolation Method: see Enzo documenation
  int method_;

};

#endif /* PROBLEM_ENZO_PROLONG_MC1_HPP */

