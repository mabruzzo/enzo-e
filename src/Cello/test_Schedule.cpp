// See LICENSE_CELLO file for license and copyright information

/// @file     test_Schedule.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Schedule class

#include "main.hpp"
#include "test.hpp"

#include "io.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ScheduleInterval");

  {
    ScheduleInterval * schedule = new ScheduleInterval;

    unit_assert (schedule != NULL);

    //--------------------------------------------------


    schedule->set_active(true);
    unit_assert(schedule->is_active());

    schedule->set_active(false);
    unit_assert(! schedule->is_active());

    schedule->set_active(true);
    unit_assert(schedule->is_active());


    unit_func("write_this_cycle");
    schedule->set_time_interval(10.0,2.0,20.0);

    unit_assert(schedule->write_this_cycle (0,  8.0) == false);
    unit_assert(schedule->write_this_cycle (0, 10.0) == true);
    unit_assert(schedule->write_this_cycle (0, 12.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 10.0) == false);
    unit_assert(schedule->write_this_cycle (0, 12.0) == true);
    unit_assert(schedule->write_this_cycle (0, 14.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 12.0) == false);
    unit_assert(schedule->write_this_cycle (0, 14.0) == true);
    unit_assert(schedule->write_this_cycle (0, 16.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 14.0) == false);
    unit_assert(schedule->write_this_cycle (0, 16.0) == true);
    unit_assert(schedule->write_this_cycle (0, 18.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 16.0) == false);
    unit_assert(schedule->write_this_cycle (0, 18.0) == true);
    unit_assert(schedule->write_this_cycle (0, 20.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 18.0) == false);
    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 22.0) == false);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0, 22.0) == false);

    unit_func("update_timestep");

    unit_assert(schedule->update_timestep(11.0, 0.1) == 0.1);
    unit_assert(schedule->update_timestep(12.0, 0.1) == 0.1);

    unit_assert (unit_incomplete);

    // schedule->set_time_interval(13.3, 2.0, 19);


    unit_func("set_time_list");

    unit_assert (unit_incomplete);
    // std::vector<double> list_double;
    // list.push_back(7/3);
    // list.push_back(-18.9);
    // list.push_back(7);
    // schedule->set_time_list (list_double);


    unit_func("this_cycle ");

    unit_assert (unit_incomplete);

    unit_func("t = schedule->update_timestep");

    unit_assert (unit_incomplete);

    unit_func("set_skip_cycle ");
    unit_assert (unit_incomplete);
    unit_func("set_skip_time ");
    unit_assert (unit_incomplete);

    delete schedule;
  }

  {
    ScheduleInterval * schedule = new ScheduleInterval;

    unit_assert (schedule != NULL);

    //--------------------------------------------------

    unit_func("write_this_cycle");
    schedule->set_minimum_time_interval(10.0,2.0,20.0);

    unit_assert(schedule->write_this_cycle (0,  8.0) == false);
    unit_assert(schedule->write_this_cycle (0, 10.0) == true);
    unit_assert(schedule->write_this_cycle (0, 12.0) == true);

    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 10.0) == false);
    unit_assert(schedule->write_this_cycle (0, 12.0) == true);
    unit_assert(schedule->write_this_cycle (0, 14.0) == true);

    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 12.0) == false);
    unit_assert(schedule->write_this_cycle (0, 14.0) == true);
    unit_assert(schedule->write_this_cycle (0, 16.0) == true);

    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 14.0) == false);
    unit_assert(schedule->write_this_cycle (0, 16.0) == true);
    unit_assert(schedule->write_this_cycle (0, 18.0) == true);

    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 16.0) == false);
    unit_assert(schedule->write_this_cycle (0, 18.0) == true);
    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 18.0) == false);
    unit_assert(schedule->write_this_cycle (0, 20.0) == true);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);
    schedule->next();

    unit_assert(schedule->write_this_cycle (0, 20.0) == false);
    unit_assert(schedule->write_this_cycle (0, 21.0) == false);


    delete schedule;

    //--------------------------------------------------

    schedule = new ScheduleInterval;
    unit_assert (schedule != NULL);
    unit_func("update_timestep");

    // for schedule_type_minimum_time, update_timestep should always just
    // return the input dt
    schedule->set_minimum_time_interval(10.0,2.0,13.0);

    unit_assert(schedule->update_timestep(9.25, 0.5) == 0.5);
    unit_assert(schedule->update_timestep(9.75, 0.5) == 0.5);
    schedule->next();
    unit_assert(schedule->update_timestep(9.25, 0.5) == 0.5);
    unit_assert(schedule->update_timestep(9.75, 0.5) == 0.5);

    delete schedule;

  }


  {
    ScheduleInterval * schedule = new ScheduleInterval;

    unit_assert (schedule != NULL);

    //--------------------------------------------------


    schedule->set_active(true);
    unit_assert(schedule->is_active());

    schedule->set_active(false);
    unit_assert(! schedule->is_active());

    schedule->set_active(true);
    unit_assert(schedule->is_active());


    unit_func("write_this_cycle");
    schedule->set_seconds_interval(1.0,3.0,5.0);

    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

    int err;
    err = system("sleep 2");

    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

    err = system("sleep 2");

    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    unit_assert(schedule->write_this_cycle (0,  0.0) == true);
    schedule->next();
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);
    unit_assert(schedule->write_this_cycle (0,  0.0) == false);

  }
  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

