// See LICENSE_CELLO file for license and copyright information
/// @file     EnzoInitialShearStream.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 19 2023
/// @brief    [\ref Enzo] Implementation of EnzoInitialShearStream

#include "enzo.hpp"
#include "charm_enzo.hpp"
#include "cello.hpp"

#include <cmath>      // std::tanh, std::exp


//----------------------------------------------------------------------

EnzoInitialShearStream::EnzoInitialShearStream(int cycle, double time)
  : Initial(cycle, time)
{}


//----------------------------------------------------------------------

static void handle_bfields_(Field& field, enzo_float target_bfield_x) {
  // this currently just handles the case where the bfield is uniform and it is
  //   - 0 along all axes (compatible with pure-hydro)
  //   - the x-axis is the only non-zero component
  bool has_bfield = (field.is_field("bfield_x") ||
                     field.is_field("bfield_y") ||
                     field.is_field("bfield_z"));
  bool has_interface_bfield = (field.is_field("bfieldi_x") ||
                               field.is_field("bfieldi_y") ||
                               field.is_field("bfieldi_z"));

  if ( !(has_bfield || has_interface_bfield) ){
    ASSERT("handle_bfields_",
           "Issue encountered while initializing shear-stream. User has "
           "specified a non-zero bfield, but hasn't initialized Enzo-E in a "
           "way that can handle bfields",
           target_bfield_x == 0);
    return;
  }

  // at this moment in time, no MHD solver can handle the case where has_bfield
  // is true & has_interface_bfield is false. Thus we assume are both true
  // right now (an error will be reported down below when we attempt to load a
  // non-existant field).
  // - In the future, we may want to revisit this if we add solvers (e.g.
  //   divergence cleaning) that don't need interface bfields

  using EFlt3DView = CelloView<enzo_float,3>;

  // define lambda function that assigns to all values in a CelloView
  auto assign_val = [](EFlt3DView arr, enzo_float val) {
    for (int iz = 0; iz < arr.shape(2); iz++){
      for (int iy = 0; iy < arr.shape(1); iy++){
        for (int ix = 0; ix < arr.shape(0); ix++){
          arr(iz,iy,ix) = val;
        }
      }
    }
  };

  // assign cell-centered values:
  assign_val(field.view<enzo_float>("bfield_x"), target_bfield_x);
  assign_val(field.view<enzo_float>("bfield_y"), 0.0);
  assign_val(field.view<enzo_float>("bfield_z"), 0.0);

  // assign face-centered values:
  assign_val(field.view<enzo_float>("bfieldi_x"), target_bfield_x);
  assign_val(field.view<enzo_float>("bfieldi_y"), 0.0);
  assign_val(field.view<enzo_float>("bfieldi_z"), 0.0);

  // now, let's update the total energy (recall: it's specific total energy)
  EFlt3DView density = field.view<enzo_float>("density");
  EFlt3DView total_energy = field.view<enzo_float>("total_energy");

  for (int iz = 0; iz < density.shape(2); iz++){
    for (int iy = 0; iy < density.shape(1); iy++){
      for (int ix = 0; ix < density.shape(0); ix++){
        total_energy(iz,iy,ix) +=
          0.5 * (target_bfield_x * target_bfield_x) / density(iz,iy,ix);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoInitialShearStream::enforce_block
(Block * block, const Hierarchy * hierarchy) throw()
{

  const double gamma = enzo::fluid_props()->gamma();

  // parameters: to be read from the parameter file
  const double chi = 100;
  const double lambda_pert = 100;
  const double rho_hot = 1.0;
  const double v_shear = 4.08;
  const double eint_hot = 10.0/(rho_hot * (gamma - 1.0)); // specific thermal
                                                          // energy
  const double vel_pert = 0.4;
  const double target_bfield_x = 0.0; // bfield component along x-axis
  const double radius = 0.5;

  const double smoothing_thickness = 0.05;
  //const double pert_loc = -2.5;
  //const double pert_width = smoothing_thickness;

  Field field = block->data()->field();
  using EFlt3DView = CelloView<enzo_float,3>;
  EFlt3DView density = field.view<enzo_float>("density");
  EFlt3DView velocity_x = field.view<enzo_float>("velocity_x");
  EFlt3DView velocity_y = field.view<enzo_float>("velocity_y");
  EFlt3DView velocity_z = field.view<enzo_float>("velocity_z");
  EFlt3DView total_energy = field.view<enzo_float>("total_energy");

  EFlt3DArray internal_energy; // specific internal (i.e. thermal) energy field
  const bool dual_energy
    = enzo::fluid_props()->dual_energy_config().any_enabled();
  if (dual_energy) {
    internal_energy = field.view<enzo_float>("total_energy");
  }

  // fetch shape of grid (including ghost zones). There aren't typos
  const int mz = density.shape(0);
  const int my = density.shape(1);
  const int mx = density.shape(2);

  // calculate x,y,z values at all cell centers (including ghost zones)
  int gx,gy,gz;
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
  std::vector<double> x_vals = std::vector<double>(mx);
  std::vector<double> y_vals = std::vector<double>(my);
  std::vector<double> z_vals = std::vector<double>(mz);
  block->data()->field_cells (x_vals.data(), y_vals.data(), z_vals.data(),
                              gx, gy, gz);

  // Now, let's fill in the grid (this loop currently includes ghost zones)
  for (int iz = 0; iz<mz; iz++){
    for (int iy = 0; iy<my; iy++){
      for (int ix = 0; ix<mx; ix++){

        const double r_cyl = std::sqrt((y_vals[iy] * y_vals[iy]) +
                                       (z_vals[iz] * z_vals[iz]));

        const double local_tanh = std::tanh((r_cyl-radius)/smoothing_thickness);

        // compute local_overdensity, relative to the wind density.
        // -> currently, it's always less than chi (nominal density contrast)
        double local_overdensity = 0.5 * (chi + 1 + (chi - 1) * -local_tanh);

        double local_density = rho_hot * local_overdensity;

        // compute the velocity component along transverse axes (y & z axes)
        double local_vy = 0.0;
        double local_vz = 0.0;

        if (lambda_pert >= 0.0) { // consider momentum perturbations

          // In the original implementation, there was a comment saying that we
          // only want to apply the perturbations to the areas between the 2
          // densities (this is a faithful transcription of the logic).
          // -> NOTE: in current implementation local_overdensity is ALWAYS
          //          smaller than chi (the nominal density contrast).
          if ((local_overdensity > 1.0) && (local_overdensity < chi)) {

            // the perturbation has a gaussian profile
            double profile = std::exp(-1 * ((r_cyl-radius) * (r_cyl-radius)) /
                                      smoothing_thickness);

            // compute perturbation along longitudinal direction (along x-axis)
            double longitude_factor = 1.0;
            if (lambda_pert > 0.0) { // the amplitude of the perturbation in
                                     // the transverse momentum component
                                     // oscilates along x-axis (as a sine-wave)
              longitude_factor = std::sin(2*cello::pi*x_vals[ix]/lambda_pert);
            } else {
              ERROR("EnzoInitialCloud::enforce_block", "uncomment this later");
              //longitude_factor =
              //  std::exp(-1 * ((x_vals[ix]-pert_loc)*(x_vals[ix]-pert_loc)) /
              //           pert_width);
            }

            // get the angle at the current location
            // - (this is like the azimuthal in cylindrical coordinates, except
            //    we have replaced z-axis with x-axis)
            // - Should we be using std::atan2(y_vals[iy], z_vals[iz])?
            double theta = std::atan(y_vals[iy]/z_vals[iz]);
            local_vz = vel_pert * profile * longitude_factor * std::cos(theta);
            local_vy = vel_pert * profile * longitude_factor * std::sin(theta);
          }

          // In the original version of the code, there was an option to add
          // random noise. We would do that here
        }

        // the following assumes uniform thermal pressure!
        double local_eint = eint_hot / local_density;

        density(iz,iy,ix) = static_cast<enzo_float>(local_density);
        velocity_x(iz,iy,ix) = static_cast<enzo_float>(v_shear * -local_tanh);
        velocity_y(iz,iy,ix) = static_cast<enzo_float>(local_vy);
        velocity_z(iz,iy,ix) = static_cast<enzo_float>(local_vz);

        if (dual_energy) {
          internal_energy(iz,iy,ix) = static_cast<enzo_float>(local_eint);
        }

        total_energy(iz,iy,ix) =
          (static_cast<enzo_float>(local_eint) +
           0.5* ( (velocity_x(iz,iy,ix) * velocity_x(iz,iy,ix)) +
                  (velocity_y(iz,iy,ix) * velocity_y(iz,iy,ix)) +
                  (velocity_z(iz,iy,ix) * velocity_z(iz,iy,ix)) ) );
      }
    }
  }

  // handle the bfields (if there are any)
  handle_bfields_(field, static_cast<enzo_float>(target_bfield_x));

  // this is where we would initialize bfields
  block->initial_done();
}
