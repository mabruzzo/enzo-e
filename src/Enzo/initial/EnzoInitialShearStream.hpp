// See LICENSE_CELLO file for license and copyright information
/// @file     EnzoInitialShearStream.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 19 2023
/// @brief    [\ref Enzo] Initialization routine for a shear stream

#ifndef ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP
#define ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP

class EnzoInitialShearStream : public Initial {
  /// @class    EnzoInitialShearStream
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialize a shearing stream. Both gas phases have
  /// bulk flows (anti-)parallel to the x-axis

public: // interface

  /// Constructor
  EnzoInitialShearStream(int cycle, double time);

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialShearStream);

  /// CHARM++ migration constructor
  EnzoInitialShearStream(CkMigrateMessage *m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;
    Initial::pup(p);
  }

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

};


#endif //ENZO_ENZO_INITIAL_SHEAR_STREAM_HPP
