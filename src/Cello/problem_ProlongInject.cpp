// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongInject.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of prolongation by injections

#include "problem.hpp"

//----------------------------------------------------------------------

ProlongInject::ProlongInject() throw() : Prolong() {
  CkPrintf("%d TRACE_PROLONG ProlongInject\n", CkMyPe());
  TRACE("ProlongInject::ProlongInject");
}

//----------------------------------------------------------------------

void ProlongInject::apply(precision_type precision, void* values_f,
                          int nd3_f[3], int im3_f[3], int n3_f[3],
                          const void* values_c, int nd3_c[3], int im3_c[3],
                          int n3_c[3], bool accumulate) {
  switch (precision) {
    case precision_single:

      apply_((float*)values_f, nd3_f, im3_f, n3_f, (const float*)values_c,
             nd3_c, im3_c, n3_c, accumulate);

      break;

    case precision_double:

      apply_((double*)values_f, nd3_f, im3_f, n3_f, (const double*)values_c,
             nd3_c, im3_c, n3_c, accumulate);

      break;

    default:

      ERROR1("ProlongInject::apply()", "Unknown precision %d", precision);
  }
}

//----------------------------------------------------------------------

template <class T>
void ProlongInject::apply_(T* values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
                           const T* values_c, int nd3_c[3], int im3_c[3],
                           int n3_c[3], bool accumulate) {
  int rank = (nd3_f[2] > 1) ? 3 : ((nd3_f[1] > 1) ? 2 : 1);

  for (int i = 0; i < rank; i++) {
    const char* xyz = "xyz";
    ASSERT3("ProlongInject::apply_",
            "fine array %c-axis %d must be 2 x coarse axis %d", xyz[i], n3_c[i],
            n3_f[i], n3_f[i] == n3_c[i] * 2);
  }

  if (n3_f[1] == 1) {
    const int ix0_f = im3_f[0];
    const int ix0_c = im3_c[0];
    const int nx_f = n3_f[0];
    if (!accumulate) {
      for (int ix_f = 0; ix_f < nx_f; ix_f++) {
        int ix_c = ix_f >> 1;

        int i_c = ix0_c + ix_c;
        int i_f = ix0_f + ix_f;

        values_f[i_f] = values_c[i_c];
      }
    } else {  // accumulate
      for (int ix_f = 0; ix_f < nx_f; ix_f++) {
        int ix_c = ix_f >> 1;

        int i_c = ix0_c + ix_c;
        int i_f = ix0_f + ix_f;

        values_f[i_f] += values_c[i_c];
      }
    }

  } else if (n3_f[2] == 1) {
    const int ix0_f = im3_f[0];
    const int iy0_f = im3_f[1];
    const int ix0_c = im3_c[0];
    const int iy0_c = im3_c[1];
    const int nx_f = n3_f[0];
    const int ny_f = n3_f[1];
    const int mx_c = nd3_c[0];
    const int mx_f = nd3_f[0];

    if (!accumulate) {
      for (int iy_f = 0; iy_f < ny_f; iy_f++) {
        int iy_c = iy_f >> 1;

        for (int ix_f = 0; ix_f < nx_f; ix_f++) {
          int ix_c = ix_f >> 1;

          int i_c = (ix0_c + ix_c) + mx_c * ((iy0_c + iy_c));
          int i_f = (ix0_f + ix_f) + mx_f * ((iy0_f + iy_f));

          values_f[i_f] = values_c[i_c];
        }
      }
    } else {  // accumulate
      for (int iy_f = 0; iy_f < ny_f; iy_f++) {
        int iy_c = iy_f >> 1;

        for (int ix_f = 0; ix_f < nx_f; ix_f++) {
          int ix_c = ix_f >> 1;

          int i_c = (ix0_c + ix_c) + mx_c * ((iy0_c + iy_c));
          int i_f = (ix0_f + ix_f) + mx_f * ((iy0_f + iy_f));

          values_f[i_f] += values_c[i_c];
        }
      }
    }

  } else {
    const int ix0_f = im3_f[0];
    const int iy0_f = im3_f[1];
    const int iz0_f = im3_f[2];
    const int ix0_c = im3_c[0];
    const int iy0_c = im3_c[1];
    const int iz0_c = im3_c[2];
    const int nx_f = n3_f[0];
    const int ny_f = n3_f[1];
    const int nz_f = n3_f[2];
    const int mx_c = nd3_c[0];
    const int mx_f = nd3_f[0];
    const int my_c = nd3_c[1];
    const int my_f = nd3_f[1];

    if (!accumulate) {
      for (int iz_f = 0; iz_f < nz_f; iz_f++) {
        int iz_c = iz_f >> 1;

        for (int iy_f = 0; iy_f < ny_f; iy_f++) {
          int iy_c = iy_f >> 1;

          for (int ix_f = 0; ix_f < nx_f; ix_f++) {
            int ix_c = ix_f >> 1;

            int i_c = (ix0_c + ix_c) +
                      mx_c * ((iy0_c + iy_c) + my_c * (iz0_c + iz_c));
            int i_f = (ix0_f + ix_f) +
                      mx_f * ((iy0_f + iy_f) + my_f * (iz0_f + iz_f));

            values_f[i_f] = values_c[i_c];
          }
        }
      }
    } else {  // accumulate
      for (int iz_f = 0; iz_f < nz_f; iz_f++) {
        int iz_c = iz_f >> 1;

        for (int iy_f = 0; iy_f < ny_f; iy_f++) {
          int iy_c = iy_f >> 1;

          for (int ix_f = 0; ix_f < nx_f; ix_f++) {
            int ix_c = ix_f >> 1;

            int i_c = (ix0_c + ix_c) +
                      mx_c * ((iy0_c + iy_c) + my_c * (iz0_c + iz_c));
            int i_f = (ix0_f + ix_f) +
                      mx_f * ((iy0_f + iy_f) + my_f * (iz0_f + iz_f));

            values_f[i_f] += values_c[i_c];
          }
        }
      }
    }
  }
}

//======================================================================
