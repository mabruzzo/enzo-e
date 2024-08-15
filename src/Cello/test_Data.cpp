// See LICENSE_CELLO file for license and copyright information

/// @file     test_Data.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar  1 20:42:57 UTC 2011
/// @brief    Test program for the Data class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN {
  PARALLEL_INIT;

  unit_init(0, 1);

  unit_class("Data");

  Data* data = new Data(3, 4, 5, 1, -1.0, 2.0, -2.0, 4.0, -3.0, 6.0,
                        new FieldDescr, new ParticleDescr);

  unit_func("Data");
  unit_assert(data != NULL);

  //----------------------------------------------------------------------

  double lower[3];
  double upper[3];

  data->lower(lower + 0, lower + 1, lower + 2);
  data->upper(upper + 0, upper + 1, upper + 2);

  unit_func("lower");
  unit_assert(lower[0] == -1.0);
  unit_assert(lower[1] == -2.0);
  unit_assert(lower[2] == -3.0);

  unit_func("upper");
  unit_assert(upper[0] == 2.0);
  unit_assert(upper[1] == 4.0);
  unit_assert(upper[2] == 6.0);

  //----------------------------------------------------------------------

  unit_func("write");

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
