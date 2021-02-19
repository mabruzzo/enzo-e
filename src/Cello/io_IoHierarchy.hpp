// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoHierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoHierarchy class

#ifndef IO_IO_HIERARCHY_HPP
#define IO_IO_HIERARCHY_HPP

class Hierarchy;

class IoHierarchy : public Io {

  /// @class    IoHierarchy
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between Hierarchy and Output classes

public: // interface

  /// Constructor
  IoHierarchy(const Hierarchy * hierarchy) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoHierarchy);

  /// CHARM++ migration constructor
  IoHierarchy(CkMigrateMessage *m) : Io(m) {}

  /// Destructor
  virtual ~IoHierarchy () throw()
  {}


  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change

    Io::pup(p);

    WARNING ("IoHierarchy::pup","skipping hierarchy_");
    //    if (p.isUnpacking()) hierarchy_ = new Hierarchy;
    //    p | *hierarchy_;
  }

#include "_io_Io_common.hpp"

  
private: // functions

  Hierarchy * hierarchy_;

private: // attributes


};

#endif /* IO_IO_HIERARCHY_HPP */

