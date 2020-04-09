// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of EnzoRiemann

#include <string>
#include <algorithm>
#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoRiemann* EnzoRiemann::construct_riemann
(std::vector<std::string> integrable_groups,
 std::vector<std::string> passive_groups, std::string solver,
 bool fallback)
{
  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(solver.size(), ' ');
  std::transform(solver.begin(), solver.end(), formatted.begin(),
		 ::tolower);
  EnzoRiemann* out;

  // Eventually we may want to check for non-MHD Riemann solvers
  if (formatted == "hll"){
    out = new EnzoRiemannHLLMHD(integrable_groups, passive_groups);
  } else if (formatted == "hlle"){
    out = new EnzoRiemannHLLEMHD(integrable_groups, passive_groups);
  } else if (formatted == "hllc"){
    out = new EnzoRiemannHLLC(integrable_groups, passive_groups);
  } else if (formatted == "hlld"){
    // could possibly check that MHD fields are included
    out = new EnzoRiemannHLLD(integrable_groups, passive_groups);
  } else {
    ERROR("EnzoRiemann::construct_riemann",
	  "The only known solvers are HLL, HLLE, HLLC, & HLLD");
    out = NULL;  // Deals with compiler warning
  }

  // May need to allow more customization of fallback solvers in the future
  if (fallback){
    ASSERT("EnzoRiemann::construct_riemann",
	   "Can't have fallback solver for the HLL & HLLE solvers",
	   formatted != "hll" && formatted != "hlle");
    out->fallback_solver_ = EnzoRiemann::construct_riemann(integrable_groups,
							   passive_groups,
							   "hlle", false);
  }

  return out;
}
