// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Initial.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar 16 09:53:31 PDT 2011
/// @brief    Implementation of the Problem class

#include "cello.hpp"
#include "main.hpp"
#include "simulation.hpp"
#include "data.hpp"
#include "problem.hpp"

//----------------------------------------------------------------------

void Initial::pup(PUP::er& p) {
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUP::able::pup(p);

  p | cycle_;
  p | time_;
}

//----------------------------------------------------------------------

void Initial::enforce_block(Block* block, const Hierarchy* hierarchy) throw() {}
