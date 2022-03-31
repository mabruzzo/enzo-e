// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodSummaryReport.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com) 
/// @date     Fri Jan 14 2022
/// @brief    [\ref Enzo] Implementation of EnzoMethodSummaryReport

#include "cello.hpp"
#include "enzo.hpp"
#include <limits> // numeric_limits

//----------------------------------------------------------------------

// these are specifically used if WRITE_TO_FILE mode

#define WRITE_TO_FILE true

#include <cstdint>
#include <cstdlib>
#include <sys/wait.h> // WEXITSTATUS

#define MAX_SUMMARY_REPORT_PER_CYCLE 3
#define N_DENS_BINS 2
#define N_EINT_BINS 50

static std::map<std::string, std::queue<std::vector<double>>> global_pending_reductions[MAX_SUMMARY_REPORT_PER_CYCLE];

static int cycle_during_calc[MAX_SUMMARY_REPORT_PER_CYCLE] = {-1, -1, -1};
static double time_during_calc[MAX_SUMMARY_REPORT_PER_CYCLE] = {-1., -1., -1.};

static std::FILE* fp_array[MAX_SUMMARY_REPORT_PER_CYCLE] =
  {nullptr, nullptr, nullptr};
static std::string* f_basename_array[MAX_SUMMARY_REPORT_PER_CYCLE] =
  {nullptr, nullptr, nullptr};

struct FileManFacade{
  // this is supposed to be light-weight and reentrant

public:
  FileManFacade(int summary_report_index)
  {
    array_index_ = summary_report_index;
    open_file_ = (fp_array[array_index_] != nullptr);
  }

  void write_entry(std::int32_t cycle, double time, double dt,
		   std::uint32_t n_dens_bins, std::uint32_t n_eint_bins,
		   std::queue<std::vector<double>>& val_q)
  {
    if (n_dens_bins != N_DENS_BINS){
      ERROR("FileManFacade::write_entry",
	    "Unexpected number of density bins");
    } else if (n_eint_bins != N_EINT_BINS){
      ERROR("FileManFacade::write_entry",
	    "Unexpected number of specific internal energy bins");
    } else if (val_q.size() != (N_DENS_BINS + N_EINT_BINS)){
      ERROR("FileManFacade::write_entry", "got the wrong number of bins");
    }

    if (dt < 0){
      ERROR("FileManFacade::write_entry",
	    "dt probably shouldn't be negative");
    }

    // open a file if one does not exist
    if (!open_file_){ open_file(time); }

    std::FILE* fp = fp_array[array_index_];

    // now actually write the entry
    {
      std::int64_t cycle_val = (std::int64_t) cycle;
      std::fwrite(&cycle_val, sizeof(std::int64_t), 1, fp);
      std::fwrite(&time, sizeof(double), 1, fp);
      std::fwrite(&dt, sizeof(double), 1, fp);
    }

    while (val_q.size() != 0){
      std::vector<double> cur_bin_vals = val_q.front();
      val_q.pop();

      std::fwrite(cur_bin_vals.data(), sizeof(double), cur_bin_vals.size(),
		  fp);
    }
  }

  void close_and_move_file(){
    // close the file
    if (!open_file_){
      ERROR("FileManFacade::close_and_move_file",
	    "can't close an already closed file");
    }
    if (std::fclose(fp_array[array_index_]) != 0){
      ERROR("FileManFacade::close_and_move_file",
	    "there was a problem closing the file");
    }
    fp_array[array_index_] = nullptr;

    // move the file
    std::string basename = *(f_basename_array[array_index_]);
    std::string cur_name = tmp_fname(basename);
    std::string target_name = "./reports/" + basename;

    // we have problems with std::rename if moving between file systems...
    std::string command = "mv " + cur_name + " " + target_name;
    if (WEXITSTATUS(std::system(command.c_str())) != 0){
      ERROR("FileManFacade::close_and_move_file",
	    "there was a problem moving the file");
    }
    delete f_basename_array[array_index_];
    f_basename_array[array_index_] = nullptr;
  }

private:
  void open_file(double time){
    if (f_basename_array[array_index_] != nullptr){
      ERROR("FileManFacade::open_file",
	    "trying to make a new file. But, f_basename_array[array_index_] "
	    "is not a nullptr");
    }
    const char* fmt = "report-t%07.4f-PEInd%d-ReportInd%d";
    int sz = std::snprintf(nullptr, 0, fmt, time, CkMyPe(), array_index_);
    std::vector<char> buf(sz + 1);
    std::snprintf(buf.data(), buf.size(), fmt, time, CkMyPe(), array_index_);

    f_basename_array[array_index_] = new std::string(buf.data());
    std::string new_fname = tmp_fname(*f_basename_array[array_index_]);

    std::FILE* fp = std::fopen(new_fname.c_str(), "wb");
    if (fp == nullptr){
      ERROR("FileManFacade::open_file",
	    "there were issues openning the file");
    } else {
      fp_array[array_index_] = fp;
    }

    // write the header
    const char name[8] = {'r',  'e',  'p',  'o',  'r',  't', '\0', '\0'};
    const char vers[8] = {'0',  '.',  '1', '\0', '\0', '\0', '\0', '\0'};
    std::int32_t n_dens_bins = N_DENS_BINS;
    std::int32_t n_eint_bins = N_EINT_BINS;

    std::fwrite(name, sizeof(char), 8, fp);
    std::fwrite(vers, sizeof(char), 8, fp);
    std::fwrite(&n_dens_bins, sizeof(std::int32_t), 1, fp);
    std::fwrite(&n_eint_bins, sizeof(std::int32_t), 1, fp);

    open_file_ = true;
  }

  std::string tmp_fname(const std::string& basename) const{
    return tmp_dir() + "/" + basename;
  }

  std::string tmp_dir() const { return "/tmp"; }

private:

  bool open_file_;
  int array_index_;

};

//----------------------------------------------------------------------

namespace{

  template<class T>
  bool map_contains_key_(const std::map<std::string,T>& m,
                         const std::string& key)
  { return m.find(key) != m.cend(); }


  // return true when a reduction is launch.
  bool launch_next_reduction_
    (Block * block,
     std::map<std::string, std::queue<std::vector<double>>>& pending_reductions)
  {
    std::string name = block->name();
    if (!map_contains_key_(pending_reductions, name)) {
      ERROR1("launch_next_reduction_",
             "pending_reductions is missing an entry for the block, %s.",
             name.c_str());
    }

    if (pending_reductions[name].empty()){
      pending_reductions.erase(name); // erase entry for name
      return false;
    } else {
      std::queue<std::vector<double>>& reduce_q = pending_reductions[name];

      // for safety, make a copy of the first element of reduce_q.
      // Then delete the first element from reduce_q
      std::vector<double> outdoubles = reduce_q.front();
      reduce_q.pop();

      EnzoBlock * enzo_block = enzo::block(block);
      CkCallback callback (CkIndex_EnzoBlock::p_method_summary_report_end(NULL),
                           enzo_block->proxy_array());
      // contribute has a special interface for passing vectors:
      enzo_block->contribute(outdoubles, CkReduction::sum_double, callback);
      return true;
    }
  }

}

//----------------------------------------------------------------------

EnzoMethodSummaryReport::EnzoMethodSummaryReport(double density_cloud,
                                                 double density_wind,
                                                 double eint_wind,
                                                 bool dual_energy,
                                                 int total_num_summaries,
                                                 int summary_report_index)
  : Method(),
    chi_(0.0),
    density_cloud_(density_cloud),
    eint_wind_(eint_wind),
    precomputed_eint_(dual_energy),
    dens_selectors_(), // we'll initialize this later when we need them
    eint_selectors_(), // we'll initialize this later when we need them
    total_num_summaries_(total_num_summaries),
    summary_report_index_(summary_report_index),
    pending_reductions_()
{
  ASSERT("EnzoMethodSummaryReport::EnzoMethodSummaryReport",
         "density_cloud, density_wind, and eint_wind must be nonzero",
         (density_cloud > density_wind) & (density_wind > 0) &
         (eint_wind > 0));
  chi_ = density_cloud / density_wind;

  ASSERT("EnzoMethodSummaryReport::EnzoMethodSummaryReport",
         "problems with total_num_summaries_/total_num_summaries",
         (total_num_summaries_ > summary_report_index_) &
         (summary_report_index_ >= 0));

  ASSERT("EnzoMethodSummaryReport::EnzoMethodSummaryReport",
	 "there are too many summary reports",
	 total_num_summaries_ <= MAX_SUMMARY_REPORT_PER_CYCLE);
}

//----------------------------------------------------------------------

namespace{
  

  template<typename T>
  struct summary_stat{
    double cell_count; // intentionally not using T to maintain precision
    T mass;
    T momentum_x;
    T energy;
    T cloud_dye_mass;

    summary_stat()
      : cell_count(0.),
        mass(0.),
        momentum_x(0.),
        energy(0.),
        cloud_dye_mass(0.)
        
    {}

    summary_stat(double cell_count, T mass, T momentum_x, T energy,
                 T cloud_dye_mass)
      : cell_count(cell_count),
        mass(mass),
        momentum_x(momentum_x),
        energy(energy),
        cloud_dye_mass(cloud_dye_mass)
    {}

    static constexpr std::size_t n_members() { return 5; }
  };



  /// actually compute the summary statistics
  ///
  /// this assumes that that ghost zones have already been clipped from each
  /// array (and that all arrays have the same shape)
  ///
  /// based on our problem size and type, we can get away with representing
  /// integer counts of cells as doubles (without any loss of precision)
  inline std::vector<summary_stat<double>> summarize_vals_
    (const std::vector<SelectFunctor_>& selectors,
     EFlt3DArray selection_field, EFlt3DArray density, EFlt3DArray velocity_x,
     EFlt3DArray eint, EFlt3DArray cloud_dye_density, double cell_volume)
    noexcept
  {
    const int nz = density.shape(0);
    const int ny = density.shape(1);
    const int nx = density.shape(2);

    std::vector<summary_stat<double>> out;
    out.reserve(selectors.size());

    for (const SelectFunctor_ selector : selectors){
      summary_stat<enzo_float> accumulator;

      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {

            bool is_selected = selector.in_bounds(selection_field(iz,iy,ix));
            enzo_float dens = is_selected * density(iz,iy,ix);
            accumulator.cell_count += is_selected;
            accumulator.mass += dens;
            accumulator.momentum_x += dens * velocity_x(iz,iy,ix);
            enzo_float cur_eint = eint(iz,iy,ix);
            // we do the extra check just for the case when we needed to
            // compute the internal energy from total_energy
            accumulator.energy += dens * cur_eint * (cur_eint >= 0);

            accumulator.cloud_dye_mass
              += is_selected * cloud_dye_density(iz,iy,ix);

          }
        }
      }

      summary_stat<double> rslt = summary_stat<double>
        (accumulator.cell_count,
         static_cast<double>(cell_volume*accumulator.mass),
         static_cast<double>(cell_volume*accumulator.momentum_x),
         static_cast<double>(cell_volume*accumulator.energy),
         static_cast<double>(cell_volume*accumulator.cloud_dye_mass));
      out.push_back(rslt);
    }
    return out;
  }
}

//----------------------------------------------------------------------

std::vector<SelectFunctor_> build_dens_selectors_(double density_cloud_,
                                                  double chi_){
  std::vector<SelectFunctor_> density_selections(2);
  density_selections[0] = {density_cloud_ / 3,
                           std::numeric_limits<enzo_float>::max()};
  density_selections[1] = {density_cloud_ / std::sqrt(chi_),
                           std::numeric_limits<enzo_float>::max()};
  return density_selections;
}

//----------------------------------------------------------------------

std::vector<SelectFunctor_> build_eint_selectors_(double eint_w,
                                                  double chi){
  double eint_cl = eint_w / chi;
  double log10_chi = std::log10(chi);

  // _tmp gives the inner bin edges in units of log(eint/eint_cl)/log(chi)
  std::vector<double> _tmp;
  if (false){ // just include bin edges
    _tmp = {0.08333333333333333, 0.25,
            0.4166666666666667,  0.5833333333333334,
            0.75, 0.9166666666666666};
  } else if (false) { // include bin edges and bin_centers
                      // (i/12 for i in range(13))
    _tmp = {0.0, 0.08333333333333333,
            0.16666666666666666, 0.25,
            0.3333333333333333, 0.4166666666666667,
            0.5, 0.5833333333333334,
            0.6666666666666666, 0.75,
            0.8333333333333334, 0.9166666666666666,
            1.0};
  } else {
    _tmp = {0.0,
	    0.020833333333333332, 0.041666666666666664, 0.0625,
	    0.08333333333333333, 0.10416666666666667, 0.125,
	    0.14583333333333334, 0.16666666666666666, 0.1875,
	    0.20833333333333334, 0.22916666666666666, 0.25,
	    0.2708333333333333, 0.2916666666666667, 0.3125,
	    0.3333333333333333, 0.3541666666666667, 0.375,
	    0.3958333333333333, 0.4166666666666667, 0.4375,
	    0.4583333333333333, 0.4791666666666667, 0.5,
	    0.5208333333333334, 0.5416666666666666, 0.5625,
	    0.5833333333333334, 0.6041666666666666, 0.625,
	    0.6458333333333334, 0.6666666666666666, 0.6875,
	    0.7083333333333334, 0.7291666666666666, 0.75,
	    0.7708333333333334, 0.7916666666666666, 0.8125,
	    0.8333333333333334, 0.8541666666666666, 0.875,
	    0.8958333333333334, 0.9166666666666666, 0.9375,
	    0.9583333333333334, 0.9791666666666666, 1.0};
  }
  std::vector<double> inner_bin_edges(_tmp.size(), 0.);
  for (std::size_t i = 0; i < inner_bin_edges.size(); i++){
    inner_bin_edges[i] = std::pow(10.0, log10_chi* _tmp[i]) * eint_cl;
  }

  std::vector<SelectFunctor_> selectors;
  selectors.reserve(inner_bin_edges.size() + 1);

  // create the selector that everything smaller than inner_bin_edges[1]
  selectors.push_back({0., (enzo_float)std::nextafter(inner_bin_edges[0],0.)});
  for (std::size_t i = 0; i < (inner_bin_edges.size()-1); i++){
    selectors.push_back({(enzo_float)(inner_bin_edges[i]),
                         (enzo_float)std::nextafter(inner_bin_edges[i+1],0.)});
  }
  // we intentionally don't use next_after on this last bin
  selectors.push_back({(enzo_float)inner_bin_edges[inner_bin_edges.size()-1],
                       std::numeric_limits<enzo_float>::max()});

  return selectors;
}

//----------------------------------------------------------------------

std::size_t EnzoMethodSummaryReport::num_selectors_ () throw()
{
  if (dens_selectors_.size() == 0){
    dens_selectors_ = build_dens_selectors_(density_cloud_, chi_);
    eint_selectors_ = build_eint_selectors_(eint_wind_, chi_);
  }

  return dens_selectors_.size() + eint_selectors_.size();
}
  

//----------------------------------------------------------------------

void EnzoMethodSummaryReport::compute ( Block * block) throw()
{
  double hx,hy,hz;
  block->data()->field_cell_width(&hx,&hy,&hz);
  double cell_volume = hx*hy*hz;

  Field field = block->data()->field();
  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);
  ASSERT("EnzoMethodSummaryReport::compute",
         "We assume that ghost depth along each axis is the same",
         (gx == gy) & (gy == gz));
  EnzoFieldArrayFactory arr_factory(block, gx);

  // get the total number of selectors (this will initialize them if they
  // haven't already been initialized)
  std::size_t expected_size = num_selectors_();

  // compute all of the summaries
  std::vector<summary_stat<double>> summaries;
  if (block->is_leaf()) {
    EFlt3DArray density = arr_factory.from_name("density");
    EFlt3DArray velocity_x = arr_factory.from_name("velocity_x");

    EFlt3DArray cloud_dye_density = arr_factory.from_name("cloud_dye");

    EFlt3DArray eint;
    if (precomputed_eint_){
      eint = arr_factory.from_name("internal_energy");
    } else {
      const int nz = density.shape(0);
      const int ny = density.shape(1);
      const int nx = density.shape(2);
      eint = EFlt3DArray(nz, ny, nx);

      EFlt3DArray velocity_y = arr_factory.from_name("velocity_y");
      EFlt3DArray velocity_z = arr_factory.from_name("velocity_z");
      EFlt3DArray etot = arr_factory.from_name("velocity_z");

      for (int iz=0; iz<nz; iz++) {
        for (int iy=0; iy<ny; iy++) {
          for (int ix=0; ix<nx; ix++) {
            enzo_float ke = 0.5* ((velocity_x(iz,iy,ix)*velocity_x(iz,iy,ix))+
                                  (velocity_y(iz,iy,ix)*velocity_y(iz,iy,ix))+
                                  (velocity_z(iz,iy,ix)*velocity_z(iz,iy,ix)));
            eint(iz,iy,ix) = etot(iz,iy,ix) - ke;
          }
        }
      }
    }

    summaries = summarize_vals_(dens_selectors_, density,
                                density, velocity_x, eint, cloud_dye_density,
                                cell_volume);

    std::vector<summary_stat<double>> tmp = summarize_vals_
      (eint_selectors_, eint,
       density, velocity_x, eint, cloud_dye_density,
       cell_volume);
    summaries.insert(summaries.end(), tmp.begin(), tmp.end());

  } else {
    summaries = std::vector<summary_stat<double>>(expected_size);
  }
  ASSERT("EnzoMethodSummaryReport::compute",
         "summaries has an unexpected length",
         summaries.size() == expected_size);

  // now, pack the summaries up and store them in reduction_queue
  std::queue<std::vector<double>> reduction_queue;

  for (const auto& summary : summaries){
    std::vector<double> tmp(summary_stat<double>::n_members());

    tmp[0] = summary.cell_count;
    tmp[1] = summary.mass;
    tmp[2] = summary.momentum_x;
    tmp[3] = summary.energy;
    tmp[4] = summary.cloud_dye_mass;

    reduction_queue.push(std::move(tmp));
  }

  // now insert the reduction_queue into pending_reductions_
  std::string block_name = block->name();
  if (map_contains_key_(pending_reductions_, block_name)) {
    ERROR1("EnzoMethodSummaryReport::compute",
           "pending_reductions_ should not already contain %s",
           block_name.c_str());
  }
  

  if (WRITE_TO_FILE){
    time_during_calc[summary_report_index_] = block->time();
    cycle_during_calc[summary_report_index_] = block->cycle();
    global_pending_reductions[summary_report_index_].emplace
      (std::move(block_name), std::move(reduction_queue));
    block->compute_done();
  } else {
    pending_reductions_.emplace(std::move(block_name),
				std::move(reduction_queue));

    // finally, launch the reduction
    bool successful_launch = launch_next_reduction_(block, pending_reductions_);

    ASSERT("EnzoMethodSummaryReport::compute", "Something went wrong",
	   successful_launch == true);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_summary_report_end(CkReductionMsg * msg)
{
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodSummaryReport::compute_resume
(Block * block,
 CkReductionMsg * msg) throw()
{
  double *values=(double *)msg->getData();
  const std::size_t n_stat_members = summary_stat<double>::n_members();
  std::size_t num_selections = (dens_selectors_.size() +
                                eint_selectors_.size());

  ASSERT("EnzoMethodSummaryReport::compute_resume",
         "long long int has unexpected size", // sanity check
         std::numeric_limits<long long int>::max() >
         std::pow(2.0,53.0) /* max integer losslessly held by double */);

  ASSERT("EnzoMethodSummaryReport::compute_resume",
         "the number of summary_stats seems to have changed!",
         n_stat_members == 5);

  if (block->index().is_root()) {

    std::string method_name = this->name();

    // determine cur_reduction_index
    std::string block_name = block->name();
    if (!map_contains_key_(pending_reductions_, block_name)) {
      ERROR1("EnzoMethodSummaryReport::compute",
             "pending_reductions_ should contain %s",
             block_name.c_str());
    }
    std::size_t remaining_reductions = pending_reductions_[block_name].size();
    std::size_t total_num_reductions = num_selectors_();
    std::size_t cur_reduction_index = // 0-index based
      total_num_reductions - (remaining_reductions + 1);

    long long int cell_count = static_cast<long long int>(values[0]);
    double mass = values[1];
    double momentum_x = values[2];
    double energy = values[3];
    double cloud_dye_mass = values[4];

    cello::monitor()->print
      ("Method",
       "%s-%d/%d- %d: %12lld, %+20.16e, %+20.16e, %+20.16e, %+20.16e",
       method_name.c_str(), summary_report_index_, total_num_summaries_,
       cur_reduction_index,
       cell_count, mass, momentum_x, energy, cloud_dye_mass);
  }

  delete msg;

  bool successful_launch = launch_next_reduction_(block, pending_reductions_);
  if (!successful_launch){
    block->compute_done();
  }
}

//----------------------------------------------------------------------

double EnzoMethodSummaryReport::timestep ( Block * block ) const throw()
{
  double out = std::numeric_limits<double>::max();

  const int in = summary_report_index_;
  if (WRITE_TO_FILE & (global_pending_reductions[in].size() != 0)){
    ASSERT("EnzoMethodSummaryReport::timestep",
           ("global_pending_reductions[summary_report_index_] must contain 0 "
	    "or 1 entries"),
	   global_pending_reductions[in].size() == 1);

    FileManFacade file_man(summary_report_index_);

    // write the output file
    std::queue<std::vector<double>> val_q;
    for (const auto& kv : global_pending_reductions[in]) {
      val_q = kv.second;
    }
    global_pending_reductions[in].clear();

    file_man.write_entry(cycle_during_calc[in], time_during_calc[in],
			 block->time() - time_during_calc[in],
			 dens_selectors_.size(), eint_selectors_.size(),
                         val_q);

    // check if we should close the file...
    double last_time_times_8 = 8*time_during_calc[in];
    int casted_last_time_times_8 = (int)last_time_times_8;
    if (last_time_times_8 == (double)casted_last_time_times_8){
      file_man.close_and_move_file();
    }
  }



  if (summary_report_index_ != 0){
    // this is necessary to avoid some problems with using multiple instances
    // in a single cycle
    return out;
  }

  // this is a somewhat hacky workaround. We would be better off updating
  // the control_stopping code
  if ( (passive_schedule_ != nullptr) &&
       (passive_schedule_->type() == schedule_type_time) ){
    ASSERT("EnzoMethodSummaryReport::timestep",
           "courant must be 1", courant_ == 1.);
    ASSERT("EnzoMethodSummaryReport::timestep",
           "courant_global must be 1",
           Method::courant_global == 1.);

    // we need to advance the schedule to the current time and cycle
    // it's ok to do this multiple times in a row as long as the arguments don't
    // change between calls
    passive_schedule_->is_scheduled(block->cycle(), block->time());

    out = passive_schedule_->time_next() - block->time();
    // SANITY CHECK
    ASSERT("MethodFrameTransform::timestep", "The timestep can't be negative.",
           out >= 0.);
  }

  return out;
}
