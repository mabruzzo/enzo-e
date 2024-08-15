// See LICENSE_CELLO file for license and copyright information

/// @file     problem_RestrictLinear.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    [\ref Problem] Declaration of the RestrictLinear class

#ifndef PROBLEM_RESTRICT_LINEAR_HPP
#define PROBLEM_RESTRICT_LINEAR_HPP

class RestrictLinear : public Restrict

{
  /// @class    RestrictLinear
  /// @ingroup  Problem
  /// @brief    [\ref Problem]

public:  // interface
  /// Constructor
  RestrictLinear() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(RestrictLinear);

  /// CHARM++ migration constructor
  RestrictLinear(CkMigrateMessage* m) : Restrict(m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p) {
    TRACEPUP;
    Restrict::pup(p);
  }

public:  // virtual functions
  /// Restrict coarse Field values to the child block (icx,icy,icz)

  virtual int apply(precision_type precision, void* values_c, int nd3_c[3],
                    int im3_c[3], int n3_c[3], const void* values_f,
                    int nd3_f[3], int im3_f[3], int n3_f[3],
                    bool accumulate = false);

  /// Return the name identifying the restrict operator
  virtual std::string name() const { return "linear"; }

private:  // functions
  template <class T>
  int apply_(T* values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
             const T* values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
             bool accumulate = false);

private:  // attributes
  // NOTE: change pup() function whenever attributes change
};

#endif /* PROBLEM_RESTRICT_LINEAR_HPP */
