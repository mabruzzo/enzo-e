// See LICENSE_CELLO file for license and copyright information

/// @file     problem_PassiveSchedule.cpp 
/// @author   Matthew Abruzzo (matthewabruzzo@gmail) 
/// @date     2020-06-22
/// @brief    [\ref Problem] Implementation of the PassiveSchedule component

#include "problem.hpp"

//----------------------------------------------------------------------

PassiveSchedule::PassiveSchedule(Schedule* schedule,
                                 double initial_time) throw()
  : PassiveSchedule()
{
  ASSERT("PassiveSchedule::PassiveSchedule",
         "the schedule pointer cannot be NULL", schedule != nullptr);
  ASSERT1("PassiveSchedule::PassiveSchedule",
          "the initial_time cannot be %.15e", PassiveSchedule::unset_time,
          initial_time != PassiveSchedule::unset_time);

  schedule_ = schedule;
  if (schedule->type() == schedule_type_minimum_time){
    int i;
    // value of limit was selected because it takes ~ 1 second, if the start
    // time is too small
    int limit = 100000000; 
    for (i=0; i < limit; i++){
      if (initial_time <= schedule->time_next()) { break; }
      schedule->next();
    }
    ASSERT2("PassiveSchedule::PassiveSchedule",
            ("Too many events (over %d) are scheduled before the "
             "initial time, %e"), limit, initial_time, (i < limit) );
    ASSERT1("PassiveSchedule::PassiveSchedule",
            "All events are scheduled for before the initial time, %e",
            initial_time,
            schedule->write_this_cycle(0, schedule->time_next()) );
  } else if (schedule->type() != schedule_type_cycle) {
    // We could potentially allow the schedule_type_time to happen, but
    // schedule_type_seconds is disallowed, unless we allow for explicit
    // synchronization
    ERROR("PassiveSchedule::PassiveSchedule",
          "The scheduling type must be either schedule_type_cycle or "
          "schedule_type_minimum_time.");
  }
}

//----------------------------------------------------------------------

PassiveSchedule::~PassiveSchedule() {
  if (schedule_) { delete schedule_; }
}

//----------------------------------------------------------------------

void PassiveSchedule::pup (PUP::er &p){
  // NOTE: update this function whenever attributes change
  p | schedule_;
  p | scheduled_event_recent_cycle_;
  p | recent_cycle_;
  p | recent_time_;
}

//----------------------------------------------------------------------

bool PassiveSchedule::is_scheduled(int cycle, double time) throw()
{
  ASSERT("PassiveSchedule::is_scheduled",
         "The wrapped schedule pointer can't be NULL", schedule_ != nullptr);
  if ((cycle != recent_cycle_) || (time != recent_time_)){
    ASSERT("PassiveSchedule::advance", "cycle and time must increase",
           (cycle > recent_cycle_) && (time>recent_time_));
    recent_cycle_ = cycle;
    recent_time_  = time;
    scheduled_event_recent_cycle_ = schedule_->write_this_cycle(cycle, time);
    if (scheduled_event_recent_cycle_){ schedule_->next(); }
  }
  return scheduled_event_recent_cycle_;
}
