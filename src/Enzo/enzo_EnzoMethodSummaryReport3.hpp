// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodSummaryReport3.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Fri Jan 28 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodSummaryReport3
///           Application specific Method to print useful diagnostic information

class EnzoMethodSummaryReport3 : public EnzoMethodSummaryReport {

  /// @class    EnzoSummaryReport
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Reports summary-statistics problems about the sim

public: // interface

  /// Create a new EnzoMethodSummaryReport object
  EnzoMethodSummaryReport3(double density_cloud, double density_wind,
                           double eint_wind, bool dual_energy,
                           int total_num_summaries, int summary_report_index)
    : EnzoMethodSummaryReport(density_cloud, density_wind, eint_wind,
                              dual_energy, total_num_summaries,
                              summary_report_index)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodSummaryReport3);

  /// Charm++ PUP::able migration constructor
  EnzoMethodSummaryReport3 (CkMigrateMessage *m)
    : EnzoMethodSummaryReport(m)
  { }

  void pup (PUP::er &p){ EnzoMethodSummaryReport::pup(p); }

  virtual void compute(Block * block) throw()
  { EnzoMethodSummaryReport::compute(block); }

  virtual std::string name () throw () 
  { return "summary_report3"; }

  virtual void compute_resume(Block * block,
                              CkReductionMsg * msg) throw()
  { EnzoMethodSummaryReport::compute_resume(block, msg); }


  virtual double timestep ( Block * block) const throw()
  { return EnzoMethodSummaryReport::timestep(block); }
};
