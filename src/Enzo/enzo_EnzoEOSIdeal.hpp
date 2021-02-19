// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the equation of state for an ideal
/// adiabatic gas

#ifndef ENZO_ENZO_EOS_IDEAL_HPP
#define ENZO_ENZO_EOS_IDEAL_HPP

class EnzoEOSIdeal : public EnzoEquationOfState
{

  /// @class    EnzoEOSIdeal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates equation of state for ideal gas
  
public: // interface
  
  /// Create a new EnzoEOSIdeal object
  EnzoEOSIdeal(double gamma, double density_floor, double pressure_floor,
               double scalar_density_floor, bool dual_energy_formalism,
	       double dual_energy_formalism_eta) throw()
    : EnzoEquationOfState(),
      gamma_(gamma),
      density_floor_(density_floor),
      pressure_floor_(pressure_floor),
      scalar_density_floor_(scalar_density_floor),
      dual_energy_formalism_(dual_energy_formalism),
      dual_energy_formalism_eta_(dual_energy_formalism_eta)
  {
    
    ASSERT("EnzoEOSIdeal", "The scalar_density_floor must not be negative.",
	   scalar_density_floor >= 0.);
    if (scalar_density_floor > 0.){
      ASSERT("EnzoEOSIdeal",
             ("scalar_density_floor is smaller than the minimum value that "
              "can be represented enzo_float"),
             scalar_density_floor >= std::numeric_limits<enzo_float>::min());
    }
  }

  /// Delete EnzoEOSIdeal object
  ~EnzoEOSIdeal()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoEOSIdeal);

  /// CHARM++ migration constructor for PUP::able
  EnzoEOSIdeal (CkMigrateMessage *m)
    : EnzoEquationOfState(m),
      gamma_(0.),
      density_floor_(0.),
      pressure_floor_(0.),
      dual_energy_formalism_(false),
      dual_energy_formalism_eta_(0.)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void reconstructable_from_integrable
  (EnzoEFltArrayMap &integrable, EnzoEFltArrayMap &reconstructable,
   EnzoEFltArrayMap &conserved_passive_map, int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const;

  void integrable_from_reconstructable
  (EnzoEFltArrayMap &reconstructable, EnzoEFltArrayMap &integrable,
   int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const;

  void pressure_from_integrable
  (EnzoEFltArrayMap &integrable_map, const EFlt3DArray &pressure,
   EnzoEFltArrayMap &conserved_passive_map, int stale_depth) const;

  void pressure_from_reconstructable(EnzoEFltArrayMap &reconstructable,
                                     EFlt3DArray &pressure,
                                     int stale_depth) const;

  inline enzo_float get_density_floor() const { return density_floor_; }

  enzo_float get_pressure_floor() const { return pressure_floor_; }

  enzo_float get_scalar_density_floor() const { return scalar_density_floor_; }

  void apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integrable_map,
                                      int stale_depth) const;

  bool is_barotropic() const { return false; }

  enzo_float get_gamma() const { return gamma_;}

  enzo_float get_isothermal_sound_speed() const { return 0;}

  // In the future, this won't be hardcoded to false
  bool uses_dual_energy_formalism() const { return dual_energy_formalism_; };

protected: // attributes
  enzo_float gamma_; // adiabatic index
  enzo_float density_floor_;
  enzo_float pressure_floor_;
  enzo_float scalar_density_floor_;
  bool dual_energy_formalism_;
  double dual_energy_formalism_eta_;
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
