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
	       bool dual_energy_formalism,
	       double dual_energy_formalism_eta) throw()
    : EnzoEquationOfState(),
      gamma_(gamma),
      density_floor_(density_floor),
      pressure_floor_(pressure_floor),
      dual_energy_formalism_(dual_energy_formalism),
      dual_energy_formalism_eta_(dual_energy_formalism_eta)
  { }

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

  void reconstructable_from_integrable(Block *block,
				       Grouping &integrable_group,
				       Grouping &reconstrable_group,
				       Grouping &conserved_passive_group,
				       int stale_depth) const;

  void integrable_from_reconstructable
  (EnzoEFltArrayMap &reconstructable, EnzoEFltArrayMap &integrable,
   int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const;

  void pressure_from_integrable(Block *block, Grouping &integrable_group,
				std::string pressure_name,
				Grouping &conserved_passive_group,
				int stale_depth) const;

  void pressure_from_reconstructable(EnzoEFltArrayMap &reconstructable,
                                     EFlt3DArray &pressure,
                                     int stale_depth) const;

  enzo_float get_density_floor() const { return density_floor_; }

  enzo_float get_pressure_floor() const { return pressure_floor_; }

  void apply_floor_to_energy_and_sync(Block *block, Grouping &integrable_group,
				      int stale_depth) const;

  bool is_barotropic() const { return false; }

  enzo_float get_gamma() const { return gamma_;}

  enzo_float get_isothermal_sound_speed() const { return 0;}

  // In the future, this won't be hardcoded to false
  bool uses_dual_energy_formalism() const { return dual_energy_formalism_; };


private:
  /// Helper function to retrieve a field array when it is possible that a
  /// field stores reconstructed data
  EFlt3DArray retrieve_field_(EnzoFieldArrayFactory &array_factory,
			      Grouping &group, std::string group_name,
			      int index, int reconstructed_axis) const;
  
  /// Copies entries of the passively advected fields included by origin_group
  /// to the corresponding entries of the fields included in destination_group
  /// reconstructed_axis = -1 means that internal field shape data can be
  /// truested. Values of 0, 1, or 2 mean that the field stores reconstructed
  /// values along the x, y, or z axis and that it is actually face-centered
  void copy_passively_advected_fields_(EnzoFieldArrayFactory &array_factory,
				       Grouping &origin_group,
				       Grouping &destination_group,
				       int reconstructed_axis = -1) const;

protected: // attributes
  enzo_float gamma_; // adiabatic index
  enzo_float density_floor_;
  enzo_float pressure_floor_;
  bool dual_energy_formalism_;
  double dual_energy_formalism_eta_;
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
