// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodSummaryReport.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Fri Jan 14 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodSummaryReport
///           Application specific Method to print useful diagnostic information

#ifndef ENZO_ENZO_METHOD_PHASE_SUMMARY_HPP
#define ENZO_ENZO_METHOD_PHASE_SUMMARY_HPP

#include <map>
#include <queue>

struct SelectFunctor_{
  enzo_float lower_bound;
  enzo_float upper_bound;

  inline bool in_bounds(enzo_float val) const noexcept {
    return (lower_bound <= val) & (val <= upper_bound);
  }
};

class EnzoMethodSummaryReport : public Method {

  /// @class    EnzoSummaryReport
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Reports summary-statistics problems about the sim

public: // interface

  /// Create a new EnzoMethodSummaryReport object
  EnzoMethodSummaryReport(double density_cloud, double density_wind,
                          double eint_wind, bool dual_energy,
                          int total_num_summaries, int summary_report_index);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodSummaryReport);

  /// Charm++ PUP::able migration constructor
  EnzoMethodSummaryReport (CkMigrateMessage *m)
    : Method (m),
      chi_(0.0),
      density_cloud_(0.0),
      eint_wind_(0.0),
      precomputed_eint_(false),
      dens_selectors_(),
      eint_selectors_(),
      total_num_summaries_(0),
      summary_report_index_(0),
      pending_reductions_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p){
    Method::pup(p);

    p | chi_;
    p | density_cloud_;
    p | eint_wind_;
    p | precomputed_eint_;
    // don't worry about pupping dens_selectors_ and eint_selectors_. We'll
    // rebuild them when we need them

    p | total_num_summaries_;
    p | summary_report_index_;

    // don't pup pending_reductions_ (It should be empty anyways)
  }

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "summary_report"; }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw();

  virtual double timestep ( Block * block) const throw();

private:

  /// returns the number of selectors. Builds the selectors if they aren't
  /// already initialized
  std::size_t num_selectors_() throw();

protected: // attributes

  double chi_;
  double density_cloud_;
  double eint_wind_;
  bool precomputed_eint_;

  std::vector<SelectFunctor_> dens_selectors_;
  std::vector<SelectFunctor_> eint_selectors_;

  /// Number of SummaryReport instances executed per cycle (per block)
  int total_num_summaries_;
  /// index of the current summary report
  int summary_report_index_;

  /// The following are used to help break up reductions into smaller message
  /// sizes
  std::map<std::string, std::queue<std::vector<double>>> pending_reductions_;
};

#endif /* ENZO_ENZO_METHOD_PHASE_SUMMARY_HPP */
