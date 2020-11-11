// See LICENSE_CELLO file for license and copyright information

/// @file     problem_PassiveSchedule.hpp 
/// @author   Matthew Abruzzo (matthewabruzzo@gmail) 
/// @date     2020-06-22
/// @brief    [\ref Problem] Declaration of the PassiveSchedule component

#ifndef PROBLEM_PASSIVE_SCHEDULE_HPP
#define PROBLEM_PASSIVE_SCHEDULE_HPP

class Schedule;

class PassiveSchedule : public PUP::able {

  /// @class    Schedule
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Schedule for events that don't actively affect
  ///           the simulation's timestep and that implicitly handles
  ///           advancement between scheduled events (and avoids explicit
  ///           synchronization).
  ///
  /// This wraps a Schedule pointer

private: // static const member

  static constexpr int unset_cycle = std::numeric_limits<int>::lowest();
  static constexpr double unset_time = std::numeric_limits<double>::lowest();

public: // functions

  /// Constructs an uninitialized version of PassiveSchedule
  PassiveSchedule() throw()
    : PUP::able(),
      schedule_(nullptr),
      scheduled_event_recent_cycle_(false),
      recent_cycle_(PassiveSchedule::unset_cycle),
      recent_time_(PassiveSchedule::unset_time)
  { }

  /// Constructs an initialized version of PassiveSchedule
  PassiveSchedule (Schedule* schedule, double initial_time) throw();

  /// can't copy without cloning wrapped pointer to schedule
  PassiveSchedule(const PassiveSchedule&) = delete;
  PassiveSchedule& operator=(const PassiveSchedule&) = delete;

  /// allow move constructor and assignment
  PassiveSchedule(PassiveSchedule&&) = default;
  PassiveSchedule& operator=(PassiveSchedule&&) = default;

  /// destructor
  ~PassiveSchedule();

  /// CHARM++ PUP::able declaration
  PUPable_decl(PassiveSchedule);

  /// CHARM++ migration constructor for PUP::able
  PassiveSchedule (CkMigrateMessage *m)
    : PUP::able(m),
      schedule_(nullptr),
      scheduled_event_recent_cycle_(false),
      recent_cycle_(PassiveSchedule::unset_cycle),
      recent_time_(PassiveSchedule::unset_time)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Return whether output is scheduled for this cycle.
  ///
  /// Consecutive calls to this method with the same arguments return the
  /// same value.
  ///
  /// The history of arguments passed to this method inform the way that the
  /// underyling schedule is advanced. In all cases, if a cycle or time is less
  /// than a value passed earlier, an unrecoverable error will be raised.
  ///
  /// @param cycle The current simulation cycle.
  /// @param time The current simulation time.
  /// @return Indicates whether there is an event scheduled during this cycle
  ///
  /// @note
  /// This accumulates a backlog of scheduled events when the simulation
  /// timestep is at least double the time interval between events
  bool is_scheduled (int cycle, double time) throw();

  /// Return the type of the underlying schedule
  int type() const throw();

  /// Return the next scheduled time
  double time_next() const throw();

private: //attributes
  /// Pointer to a schedule object
  Schedule *schedule_;

  /// Indicates whether a scheduled event occured in the recent cycle
  bool scheduled_event_recent_cycle_;

  /// Records the most recent cycle when advance was invoked
  int recent_cycle_;

  /// Records the most recent time when advance was invoked
  double recent_time_;
};

#endif /* PROBLEM_PASSIVE_SCHEDULE_HPP */
