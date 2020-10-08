// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImpl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of RiemannImpl, which is a class
/// template that can be specialized to implement various Riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL_HPP
#define ENZO_ENZO_RIEMANN_IMPL_HPP

#include <pup_stl.h>
#include <cstdint> // used to check that static methods are defined

//----------------------------------------------------------------------

struct HydroLUT {
  /// @class    HydroLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates a compile-time LUT for pure
  ///           hydrodynamics that is used for implementing Riemann Solvers

  enum vals { density=0,
	      velocity_i,
	      velocity_j,
	      velocity_k,
	      total_energy,
	      NEQ};

  // in the future, automatically calculate the following
  static const std::size_t specific_start = 1;
};

//----------------------------------------------------------------------

struct MHDLUT{
  /// @class    MHDLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates a compile-time LUT for MHD that is
  ///           that is to be used for implementing Riemann Solvers
  enum vals { density=0,
	      bfield_i,
	      bfield_j,
	      bfield_k,
	      velocity_i,
	      velocity_j,
	      velocity_k,
	      total_energy,
	      NEQ};
  // in the future, automatically calculate the following
  static const std::size_t specific_start = 4;
};

//----------------------------------------------------------------------

/// @typedef riemann_function_call_signature
/// @brief   This is the function call signature of the function call operator,
///          `operator()`, that is expected for an `ImplFunctor` used to
///          implement a Riemann Solver by specializing `EnzoRiemannImpl`.
///
/// An `ImplFunctor`'s `operator()` method should be declared as:
///
/// @code
///     lutarray<LUT> operator()
///       (const lutarray<LUT> flux_l, const lutarray<LUT> flux_r,
///        const lutarray<LUT> prim_l, const lutarray<LUT> prim_r,
///        const lutarray<LUT> cons_l, const lutarray<LUT> cons_r,
///        enzo_float pressure_l, enzo_float pressure_r, bool barotropic_eos,
///        enzo_float gamma, enzo_float isothermal_cs, enzo_float &vi_bar)
/// @endcode
///
/// This function should computes the Riemann Flux at a given cell-interface
/// for the set of actively advected quantities designated by the `LUT` type
/// publicly defined within the scope of `ImplFunctor` (that can be acessed with
/// `ImplFunctor::LUT`).
///
/// The arguments, `flux_l`, `flux_r`, `prim_l`, `prim_r`, `cons_l`, and
/// `cons_r` are constant arrays that store values associated with each
/// actively advected quantity, at the indices designated by
/// `ImplFunctor::LUT`. For each quantity, the arrays hold the its associated
/// fluxes, its inherent values, and its values when converted to conserved
/// form. Note that because the quantities are inherently specific or
/// conserved. For quantities that are inherently conserved (e.g. density) the
/// associated values in `prim_l`/`prim_r` are identical those in
/// `cons_l`/`cons_r`.
///
/// The left and right reconstructed pressure values are passed as `pressure_l`
/// and `pressure_r`. `barotropic_eos` indicates whether the fluid equation of
/// state is barotropic. If `true`, then `isothermal_cs` is expected to be
/// non-zero and if `false`, then `gamma` is expected to be positive. `vi_bar`
/// is the estimate of `velocity_i` at the cell-interface to be used in the
/// internal energy source term (for the dual energy formalism).
///
/// The function is expected to return the computed fluxes for each actively
/// advected quantity in an array with mapping set by `ImplFunctor::LUT`
template<class LUT>
using riemann_function_call_signature =
  lutarray<LUT>(*)(const lutarray<LUT> flux_l, const lutarray<LUT> flux_r,
                   const lutarray<LUT> prim_l, const lutarray<LUT> prim_r,
                   const lutarray<LUT> cons_l, const lutarray<LUT> cons_r,
                   enzo_float pressure_l, enzo_float pressure_r,
                   bool barotropic_eos, enzo_float gamma,
                   enzo_float isothermal_cs, enzo_float &vi_bar);

//----------------------------------------------------------------------

template <class ImplFunctor>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam ImplFunctor The functor used to specialize `EnzoRiemannImpl`. The
  ///     functor must provide a public member type called `LUT`, which is a
  ///     specialization of `EnzoRiemannLUT<InputLUT>`. It must also be default
  ///     constructible and provide a public `operator()` method. For details
  ///     about the method's expected signature, see the documentation for
  ///     `riemann_function_call_signature`.
  ///
  /// EnzoRiemannImpl factors out the repeated code between different
  /// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF &
  /// Roe solvers).
  ///
  /// To define a new RiemannSolver using `EnzoRiemann`:
  ///   1. Define a new `ImplFunctor` (e.g. `HLLDImpl`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoRiemannImpl` that uses the new `ImplFunctor`
  ///      (e.g. `using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;`).
  ///   3. Then add the particlular specialization to enzo.CI (e.g. add the
  ///      line: `PUPable EnzoRiemannImpl<HLLDImpl>;`)
  ///   4. Update `EnzoRiemann::construct_riemann` to construct the Riemann
  ///      Solver when the correct name is specified.
  ///   5. Update the documentation with the name of the newly available
  ///      RiemannSolver

  using LUT = typename ImplFunctor::LUT;

  static_assert(std::is_default_constructible<ImplFunctor>::value,
		"ImplFunctor is not default constructable");

  // Check whether ImplFunctor's operator() method has the expected signature
  // and raise an error message if it doesn't:
  //   - first, define a struct (has_expected_functor_sig_) to check if
  //     ImplFunctor's operator() method has the expected signature. This needs
  //     to happen in the current scope because the signature depends on LUT
  DEFINE_HAS_INSTANCE_METHOD(has_expected_functor_sig_, operator(),
                             riemann_function_call_signature<LUT>);
  static_assert(has_expected_functor_sig_<ImplFunctor>::value,
		"ImplFunctor's operator() method doesn't have the correct "
		"function signature");

  /// The following vector holds advection quantities for which passive
  /// advection fluxes can be computed when they are not directly handled by a
  /// Riemann Solver. Unlike normal passively advected scalars, quantities
  /// listed here nominally play a more direct role in hydrodynamical
  /// integration. The canonical example is "specific internal energy"
  ///
  /// Note: This is unrelated to Riemann Solver Fallback. 
  static const std::vector<std::string> PassiveFallbackAdvectionQuantities;

public: // interface

  /// Constructor
  ///
  /// @param integrable_groups A vector of integrable quantities (listed as
  ///     advected quantities in FIELD_TABLE). These are used as group names in
  ///     the Grouping objects that store field names. In effect this is used
  ///     to register the quantities operated on by the Riemann Solver
  /// @param passive_groups A vector with the names of the groups of passively
  ///     advected scalars that may be included. (If a group is listed here but
  ///     the Grouping object doesn't actually provide any fields in the group,
  ///     no problems are caused
  EnzoRiemannImpl(std::vector<std::string> integrable_groups,
                   std::vector<std::string> passive_groups);

  /// Virtual destructor
  virtual ~EnzoRiemannImpl(){ };

  /// CHARM++ PUP::able declaration
  PUPable_decl_template(EnzoRiemannImpl<ImplFunctor>);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemannImpl (CkMigrateMessage *m)
    : EnzoRiemann(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  ///TEMPORARY
  /// This is to assist with helping the transition from Groupings to maps
  const std::vector<std::string> combined_integrable_groups() const throw();

  void solve (Block *block, Grouping &priml_group, Grouping &primr_group,
	      std::string pressure_name_l, std::string pressure_name_r,
	      Grouping &flux_group, int dim, EnzoEquationOfState *eos,
	      int stale_depth, std::string interface_velocity_name) const;

protected : //methods
  
  /// Computes the fluxes for the passively advected quantites.
  void solve_passive_advection_(Block* block, Grouping &priml_group,
				Grouping &primr_group, Grouping &flux_group,
				const EFlt3DArray &density_flux, int dim,
				int stale_depth) const throw();

protected: //attributes

  /// Names of the quantities to advect
  std::vector<std::string> integrable_groups_;

  /// Names of the passively advected groups of fields (e.g. colours)
  std::vector<std::string> passive_scalar_groups_;

  /// Names of the integrable quantities that require passive fluxes
  std::vector<std::string> passive_integrable_groups_;

  /// Field categories of each quantity in passive_integrable_groups_
  std::vector<FieldCat> passive_integrable_categories_;

  // TODO: Remove the following attribute. This is merely a stop-gap solution
  //       to support the use of EnzoRiemannHLLC with VL+CT (in the future, it
  //       will be possible to run CL integrator with magnetic fields).
  /// Indicates whether the bfields should be forced to zero
  bool force_bfield_fluxes_to_zero_;
};

//----------------------------------------------------------------------

// TODO: move this to EnzoRiemann so that it isn't initialized for every
// template

template <class ImplFunctor>
const std::vector<std::string>
EnzoRiemannImpl<ImplFunctor>::PassiveFallbackAdvectionQuantities =
  {"internal_energy"};

//----------------------------------------------------------------------

template <class ImplFunctor>
EnzoRiemannImpl<ImplFunctor>::EnzoRiemannImpl
(std::vector<std::string> integrable_groups,
 std::vector<std::string> passive_groups)
  : EnzoRiemann()
{ 
  // Quick sanity check - integrable_groups must have density and velocity
  ASSERT("EnzoRiemannImpl","integrable_groups must contain \"density\"",
	 std::find(integrable_groups.begin(), integrable_groups.end(),
		   "density") != integrable_groups.end());
  ASSERT("EnzoRiemannImpl","integrable_groups must contain \"velocity\"",
	 std::find(integrable_groups.begin(), integrable_groups.end(),
		   "velocity") != integrable_groups.end());

  LUT::validate(); // sanity check - validate LUT's expected properties

  // Determine all quantities used by LUT
  std::set<std::string> lut_groups = LUT::quantity_names();

  // TODO: add special treatment of total energy for cases with barotropic EOS

  // Check that all quantities in lut_groups are in integrable_groups
  for (std::string group : lut_groups){
    ASSERT1("EnzoRiemannImpl", "\"%s\" must be a registered integrable group",
	    group.c_str(),
	    std::find(integrable_groups.begin(), integrable_groups.end(),
		      group) != integrable_groups.end());
  }

  force_bfield_fluxes_to_zero_ = false;

  // Check that all quantities in integrable_groups are in lut_groups or in
  // EnzoRiemannImpl<ImplFunctor>::PassiveFallbackAdvectionQuantities
  for (std::string group : integrable_groups){
    if (std::find(lut_groups.begin(), lut_groups.end(), group)
	!= lut_groups.end()){
      integrable_groups_.push_back(group);
    } else if (std::find(PassiveFallbackAdvectionQuantities.begin(),
			 PassiveFallbackAdvectionQuantities.end(),
			 group) != PassiveFallbackAdvectionQuantities.end()){
      passive_integrable_groups_.push_back(group);
    } else if (group == "bfield") {
      WARNING("EnzoRiemannImpl",
              "This solver isn't equipped to handle bfields. Their fluxes "
              "will all be set to 0.");
      force_bfield_fluxes_to_zero_ = true;
    } else {
      ERROR1("EnzoRiemannImpl", "can't handle the \"%s\" integrable group",
	     group.c_str());
    }
  }

  // determine the categories of each passive integrable group
  for (std::string group : passive_integrable_groups_){
    FieldCat category;
    EnzoCenteredFieldRegistry::quantity_properties(group, nullptr,
                                                   &category, nullptr);
    passive_integrable_categories_.push_back(category);
  }

  passive_scalar_groups_ = passive_groups;
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::pup (PUP::er &p)
{
  EnzoRiemann::pup(p);

  p|integrable_groups_;
  p|passive_scalar_groups_;
  p|passive_integrable_groups_;
  p|passive_integrable_categories_;
  p|force_bfield_fluxes_to_zero_;
}

//----------------------------------------------------------------------

template <class ImplFunctor>
const std::vector<std::string> EnzoRiemannImpl<ImplFunctor>::combined_integrable_groups() const throw()
{
  std::vector<std::string> combined_integrable_groups;

  auto contains = [](const std::vector<std::string> &vec,
		     const std::string &value)
  { return std::find(vec.cbegin(),vec.cend(), value) != vec.cend(); };

  for (const std::string& elem : integrable_groups_){
    if (!contains(combined_integrable_groups, elem)) {
      combined_integrable_groups.push_back(elem);
    }
  }
  for (const std::string& elem : passive_scalar_groups_){
    if (!contains(combined_integrable_groups, elem)) {
      combined_integrable_groups.push_back(elem);
    }
  }
  for (const std::string& elem : passive_integrable_groups_){
    if (!contains(combined_integrable_groups, elem)) {
      combined_integrable_groups.push_back(elem);
    }
  }

  return combined_integrable_groups;
}

//----------------------------------------------------------------------

inline void set_bfield_fluxes_to_zero_(EnzoFieldArrayFactory array_factory,
                                       Grouping &flux_group) noexcept
{
  int num_fields = flux_group.size("bfield");
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    EFlt3DArray bfield_flux = array_factory.from_grouping(flux_group, "bfield",
                                                          field_ind);
    for (int iz=0; iz<bfield_flux.shape(0); iz++) {
      for (int iy=0; iy<bfield_flux.shape(1); iy++) {
        for (int ix=0; ix<bfield_flux.shape(2); ix++) {
          bfield_flux(iz,iy,ix) = 0.;
        }
      }
    }

  }
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve
(Block *block, Grouping &priml_group, Grouping &primr_group,
 std::string pressure_name_l, std::string pressure_name_r,
 Grouping &flux_group, int dim, EnzoEquationOfState *eos, int stale_depth,
 std::string interface_velocity_name) const
{

  const bool barotropic = eos->is_barotropic();
  // if not barotropic then the following doesn't have to be reasonable
  const enzo_float isothermal_cs = eos->get_isothermal_sound_speed();
  // If barotropic, then the following doesn't have to be a reasonable value
  // (this will have to be tweaked when we introduce species that modify gamma)
  const enzo_float gamma = eos->get_gamma();

  // When barotropic equations of state are eventually introduced, all eos
  // dependencies should be moved up here
  ASSERT("EnzoRiemannImpl::solve", "currently no support for barotropic eos",
	 !barotropic);

  // Load arrays:
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  // First, load in the precomputed pressure 
  EFlt3DArray pressure_array_l, pressure_array_r;
  pressure_array_l = array_factory.assigned_center_from_name(pressure_name_l,
							     dim);
  pressure_array_r = array_factory.assigned_center_from_name(pressure_name_r,
							     dim);

  // Second, load field to store interface velocity (if applicable)
  Field field = block->data()->field();
  EFlt3DArray velocity_i_bar_array;
  const bool store_interface_vel =((interface_velocity_name != "") ||
				   (field.is_field(interface_velocity_name)));

  if (store_interface_vel) {
    velocity_i_bar_array =
      array_factory.assigned_center_from_name(interface_velocity_name, dim);
  }

  std::array<EFlt3DArray,LUT::NEQ> wl_arrays, wr_arrays, flux_arrays;
  using enzo_riemann_utils::load_array_of_fields;
  wl_arrays = load_array_of_fields<LUT>(block, priml_group, dim, stale_depth);
  wr_arrays = load_array_of_fields<LUT>(block, primr_group, dim, stale_depth);
  flux_arrays = load_array_of_fields<LUT>(block, flux_group, dim,
					  stale_depth);

  ImplFunctor func;

  // For Nearest-Neighbor, we care about interfaces starting at i+1/2. Thanks
  // to use of stale_depth, we can treat every case like this
  for (int iz=0; iz<flux_arrays[0].shape(0); iz++) {
    for (int iy=0; iy<flux_arrays[0].shape(1); iy++) {
      for (int ix=0; ix<flux_arrays[0].shape(2); ix++) {

	lutarray<LUT> wl, wr;
	// get the fluid fields
	for (std::size_t field_ind=0; field_ind<LUT::NEQ; field_ind++){
	  wl[field_ind] = wl_arrays[field_ind](iz,iy,ix);
	  wr[field_ind] = wr_arrays[field_ind](iz,iy,ix);
	}

	// get the left/right pressure
	enzo_float pressure_l = pressure_array_l(iz,iy,ix);
	enzo_float pressure_r = pressure_array_r(iz,iy,ix);

	// get the conserved quantities
	lutarray<LUT> Ul = enzo_riemann_utils::compute_conserved<LUT>(wl);
	lutarray<LUT> Ur = enzo_riemann_utils::compute_conserved<LUT>(wr);

	// compute the interface fluxes
	lutarray<LUT> Fl = enzo_riemann_utils::active_fluxes<LUT>(wl, Ul,
                                                                pressure_l);
	lutarray<LUT> Fr = enzo_riemann_utils::active_fluxes<LUT>(wr, Ur,
                                                                pressure_r);

	enzo_float interface_velocity_i;
	// Now compute the Riemann Fluxes
        lutarray<LUT> fluxes = func(Fl, Fr, wl, wr, Ul, Ur, pressure_l,
                                    pressure_r, barotropic, gamma,
                                    isothermal_cs, interface_velocity_i);

	// record the Riemann Fluxes
	for (std::size_t field_ind=0; field_ind<LUT::NEQ; field_ind++){
	  flux_arrays[field_ind](iz,iy,ix) = fluxes[field_ind];
	}

	if (store_interface_vel){
	  velocity_i_bar_array(iz,iy,ix) = interface_velocity_i;
	}
      }
    }
  }

  solve_passive_advection_(block, priml_group, primr_group, flux_group,
  			   flux_arrays[LUT::density], dim, stale_depth);

  if (LUT::has_bfields()){
    // If Dedner Fluxes are required, they might get handled here
    //   - It would probably be better to handle this separately from the
    //     Riemann Solver since we already precompute L & R conserved
    //     AND it doesn't require wavespeed information.
  } else {
    // TODO: finish refactoring the VL+CT method so it can be ran entirely
    //        without magnetic fields (so that we can remove the following):
    if (force_bfield_fluxes_to_zero_){
      set_bfield_fluxes_to_zero_(array_factory, flux_group);
    }
  }
}

//----------------------------------------------------------------------

inline void compute_unity_sum_passive_fluxes_(const enzo_float dens_flux,
					      EFlt3DArray *flux_arrays,
					      EFlt3DArray *reconstructed,
					      const int num_fields,
					      const int iz, const int iy,
					      const int ix) noexcept
{
  ERROR("EnzoRiemannImpl::solve_passive_advection_",
        "This codepath is untested");
  enzo_float sum = 0.;
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    sum += reconstructed[field_ind](iz,iy,ix);
  }
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    flux_arrays[field_ind](iz,iy,ix) = (reconstructed[field_ind](iz,iy,ix)
					* dens_flux / sum);
  }
}

//----------------------------------------------------------------------

inline void passive_advection_helper_
(std::string group_name, Grouping &priml_group, Grouping &primr_group,
 Grouping &flux_group, const EFlt3DArray &density_flux, int dim,
 EnzoFieldArrayFactory array_factory, const bool unity_sum) noexcept
{
  int num_fields = priml_group.size(group_name);
  if (num_fields == 0) {return;}
  EnzoPermutedCoordinates coord(dim);

  // This was basically transcribed from hydro_rk in Enzo:

  // load array of fields
  EFlt3DArray *wl_arrays = new EFlt3DArray[num_fields];
  EFlt3DArray *wr_arrays = new EFlt3DArray[num_fields];
  EFlt3DArray *flux_arrays = new EFlt3DArray[num_fields];
  
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    wl_arrays[field_ind] = array_factory.reconstructed_field(priml_group,
							     group_name,
							     field_ind, dim);
    wr_arrays[field_ind] = array_factory.reconstructed_field(primr_group,
							     group_name,
							     field_ind, dim);
    flux_arrays[field_ind] = array_factory.from_grouping(flux_group,
							 group_name,
							 field_ind);
  }

  for (int iz=0; iz<density_flux.shape(0); iz++) {
    for (int iy=0; iy<density_flux.shape(1); iy++) {
      for (int ix=0; ix<density_flux.shape(2); ix++) {

	enzo_float dens_flux = density_flux(iz,iy,ix);
	EFlt3DArray *reconstr = (dens_flux>0) ? wl_arrays : wr_arrays;

	if (unity_sum) {
	  compute_unity_sum_passive_fluxes_(dens_flux, flux_arrays, reconstr,
					     num_fields, iz, iy, ix);
	} else {
	  for (int field_ind=0; field_ind<num_fields; field_ind++){
	    flux_arrays[field_ind](iz,iy,ix) = (reconstr[field_ind](iz,iy,ix)
						* dens_flux);
	  }
	}

      }
    }
  }
  delete[] wl_arrays; delete[] wr_arrays; delete[] flux_arrays;
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve_passive_advection_
(Block* block, Grouping &priml_group, Grouping &primr_group,
 Grouping &flux_group, const EFlt3DArray &density_flux, int dim,
 int stale_depth)
  const throw()
{

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  for (std::string group_name : passive_scalar_groups_){
    // If we know of group_name that must always sum to 1, we should check for
    // that now! (The easiest way to facillitate that would be to separate
    // species from colours)
    bool unity_sum = false;

    // it's OK if no fields were included in a given group because the included
    // fields may be entirely determined by the methods unrelated to the hydro
    // integration method
    passive_advection_helper_(group_name, priml_group,primr_group, flux_group,
			      density_flux, dim, array_factory, unity_sum);
  }

  for (std::size_t i=0; i<passive_integrable_groups_.size(); i++){
    std::string group_name = passive_integrable_groups_[i];

    ASSERT1("EnzoRiemannImpl::solve_passive_advection_",
	    "Passive advection fluxes were to be computed for the integrable "
	    "quantity, \"%s\", but no fields were provided.",
	    group_name.c_str(),  priml_group.size(group_name) != 0);

    if (passive_integrable_categories_[i] != FieldCat::specific){
      // this is here to indicate that FieldCat::conserved should be allowed
      INCOMPLETE("passive_advection_helper_ doesn't currently support "
                 "conserved passive scalars");
    }
    passive_advection_helper_(group_name, priml_group,primr_group, flux_group,
			      density_flux, dim, array_factory, false);
  }
}

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
