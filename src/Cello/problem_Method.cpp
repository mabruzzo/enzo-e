// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-09-04
/// @brief    

#include "problem.hpp"

double Method::courant_global = 1.0;

//----------------------------------------------------------------------

Method::~Method() throw()
{
  delete passive_schedule_;
  for (size_t i=0; i<refresh_list_.size(); i++) {
    delete refresh_list_[i];
    refresh_list_[i] = 0;
  }
}

//----------------------------------------------------------------------

void Method::pup (PUP::er &p)
{ TRACEPUP;
  PUP::able::pup(p);

  bool pk = p.isPacking();
  bool up = p.isUnpacking();

  int n;
  if (pk) n=refresh_list_.size();
  p | n;
  if (up) refresh_list_.resize(n);
  for (int i=0; i<n; i++) {
    p | refresh_list_[i]; // PUP::able
  }

  p | refresh_list_;
  p | passive_schedule_; // pupable
  p | courant_;
}

//----------------------------------------------------------------------

bool Method::advance_schedule(int cycle, double time) throw(){
  if (passive_schedule_){
    // Implicitly handles advancement
    return passive_schedule_->is_scheduled(cycle,time);
  } else { // default case: scheduled every cycle
    return true;
  }
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule, double initial_time) throw()
{ 
  if (passive_schedule_) delete passive_schedule_;
  passive_schedule_ = new PassiveSchedule(schedule, initial_time);
}

//======================================================================

