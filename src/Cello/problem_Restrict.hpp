// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Restrict.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    [\ref Problem] Declaration of the Restrict class
///

#ifndef PROBLEM_RESTRICT_HPP
#define PROBLEM_RESTRICT_HPP

class Restrict : public PUP::able

{
  /// @class    Restrict
  /// @ingroup  Problem
  /// @brief    [\ref Problem]

public:  // interface
  /// Constructor
  Restrict() throw();

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Restrict);

  Restrict(CkMigrateMessage* m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p) {
    TRACEPUP;
    PUP::able::pup(p);
  }

public:  // virtual functions
  /// Restrict coarse Field values to the child block (icx,icy,icz)

  virtual int apply(precision_type precision, void* values_c, int nd3_c[3],
                    int im3_c[3], int n3_c[3], const void* values_f,
                    int nd3_f[3], int im3_f[3], int n3_f[3],
                    bool accumulate = false) = 0;

  /// Return the name identifying the restrict operator
  virtual std::string name() const = 0;

private:  // functions
private:  // attributes
  // NOTE: change pup() function whenever attributes change
};

#endif /* PROBLEM_RESTRICT_HPP */
