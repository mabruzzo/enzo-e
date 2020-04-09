// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon June 24 2019
/// @brief    [\ref Enzo] Implementation of EnzoIntegrableUpdate.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoIntegrableUpdate::EnzoIntegrableUpdate
(std::vector<std::string> integrable_groups, bool skip_B_update,
 std::vector<std::string> passive_groups) throw()
{
  integrable_groups_ = integrable_groups;
  passive_groups_ = passive_groups;

  if (skip_B_update){
    // remove bfield group from integrable_groups (if present)
    std::string target("bfield");
    integrable_groups_.erase(std::remove(integrable_groups_.begin(),
					 integrable_groups_.end(), target),
			     integrable_groups_.end());
  }
  setup_lut_();

  // sanity check:
  std::vector<std::string>::size_type num_unique_groups;
  num_unique_groups = this->combined_integrable_groups().size();
  ASSERT("EnzoIntegrableUpdate",
	 ("group names appear more than once in integrable_groups or "
	  "passive_groups"),
	 num_unique_groups == (integrable_groups_.size() +
			       passive_groups_.size()) );
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::setup_lut_(){
  EnzoCenteredFieldRegistry registry;
  lut_ = registry.prepare_advection_lut(integrable_groups_,
					conserved_start_, conserved_stop_,
					specific_start_, specific_stop_,
				        other_start_, other_stop_, nfields_);

  ASSERT("EnzoMethodMHDVlct::update_quantities_",
	 "not equipped to handle fields classified as other",
	 other_start_ == other_stop_);
}

//----------------------------------------------------------------------

const std::vector<std::string> EnzoIntegrableUpdate::combined_integrable_groups
(bool omit_flagged) const throw()
{
    auto contains = [](const std::vector<std::string> &vec,
			const std::string &value) -> bool
      { return std::find(vec.cbegin(),vec.cend(), value) != vec.cend(); };

    std::vector<std::string> out;
    for (const std::string& elem : integrable_groups_){
      if (!contains(out, elem)) { out.push_back(elem); }
    }
    for (const std::string& elem : passive_groups_){
      if (!contains(out, elem)) { out.push_back(elem); }
    }
    return out;
}

//----------------------------------------------------------------------

void get_limits_(EnzoFieldArrayFactory &af,
		 std::vector<std::string> group_names, Grouping &dUcons_group,
		 int *nx, int *ny, int *nz){
  for (const std::string& name : group_names){
    int num_fields = dUcons_group.size(name);
    if (num_fields > 0){
      EFlt3DArray arr = af.from_grouping(dUcons_group, name, 0);
      *nz = arr.shape(0);
      *ny = arr.shape(1);
      *nx = arr.shape(2);
      break;
    }
  }
}

//----------------------------------------------------------------------

void check_fallback_mask_shape_(EnzoFieldArrayFactory &af,
				std::string function_name,
				std::vector<std::string> group_names,
				Grouping &dUcons_group,
				CelloArray<bool,3> &mask)
{
  int nx, ny, nz;
  get_limits_(af, group_names,dUcons_group, &nx, &ny, &nz);
  ASSERT(function_name.c_str(),
	 "Field in dUcons_group and fallback_mask don't have same shape\n",
	 (nz == mask.shape(0) && ny == mask.shape(1) && nx == mask.shape(2)));
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::clear_dUcons_group
(Block *block, Grouping &dUcons_group, enzo_float value,
 CelloArray<bool,3> *fallback_mask_ptr) const
{
  const std::vector<std::string> group_names =combined_integrable_groups(true);
  EnzoFieldArrayFactory array_factory(block, 0);

  const bool masked = fallback_mask_ptr != NULL;

  CelloArray<bool,3> mask;
  if (masked){
    mask = *fallback_mask_ptr;
    check_fallback_mask_shape_(array_factory,
			       "EnzoIntegrableUpdate::clear_dUcons_group",
			       group_names, dUcons_group, mask);
  }
  
  for (const std::string& name : group_names){
    int num_fields = dUcons_group.size(name);

    for (int i=0; i<num_fields; i++){
      EFlt3DArray array = array_factory.from_grouping(dUcons_group, name, i);

      for (int iz=0; iz<array.shape(0); iz++) {
	for (int iy=0; iy<array.shape(1); iy++) {
	  for (int ix=0; ix<array.shape(2); ix++) {
	    if (!masked || mask(iz,iy,ix)) { array(iz,iy,ix) = value; }
	  }
	}
      }

    }
  }
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::accumulate_flux_component
(Block *block, int dim, double dt, Grouping &flux_group,
 Grouping &dUcons_group, int stale_depth,
 CelloArray<bool,3> *fallback_mask_ptr) const
{
  const std::vector<std::string> group_names =combined_integrable_groups(true);
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EnzoPermutedCoordinates coord(dim);

  CSlice full_ax(nullptr, nullptr);

  const bool masked = fallback_mask_ptr != NULL;
  CelloArray<bool,3> mask;
  if (masked){
    CSlice slc(stale_depth, -1*stale_depth);
    CelloArray<bool,3> temp_mask = fallback_mask_ptr->subarray(slc, slc, slc);
    check_fallback_mask_shape_
      (array_factory, "EnzoIntegrableUpdate::accumulate_flux_component",
       group_names, dUcons_group, temp_mask);
    mask = coord.get_subarray(temp_mask, full_ax, full_ax, CSlice(1, -1));
  }

  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx_i = dt/enzo_block->CellWidth[coord.i_axis()];

  

  for (std::string name : group_names){

    int num_fields = dUcons_group.size(name);
    ASSERT1("EnzoIntegrableUpdate::accumulate_flux_divergence",
	    ("The number of fields in the \"%s\" group is differnt between "
	     "flux_group and dUcons_group."), name.c_str(),
	    num_fields == flux_group.size(name));

    for (int field_ind = 0; field_ind < num_fields; field_ind++){
      // Since we don't have fluxes on the exterior faces along axis i, we can
      // not update dU in the first & last cell along the axis. Thus we cast
      // the calculation as:
      //   dU(k,j,i+1) -= dtdx_i * (flux(k, j, i+3/2) - flux(k, j, i+1/2))

      // define : dU_center(k,j,i) -> dU(k,j,i+1)
      EFlt3DArray dU, dU_center;
      dU = array_factory.from_grouping(dUcons_group, name, field_ind);
      dU_center = coord.get_subarray(dU, full_ax, full_ax, CSlice(1, -1));

      // define:  fl(k,j,i)        -> flux(k, j, i+1/2)
      //          fr(k,j,i)        -> flux(k, j, i+3/2)
      EFlt3DArray flux, fl, fr;
      flux = array_factory.from_grouping(flux_group, name, field_ind);
      fl = coord.get_subarray(flux, full_ax, full_ax, CSlice(0, -1));
      fr = coord.get_subarray(flux, full_ax, full_ax, CSlice(1, nullptr));

      auto inner_loop = [dtdx_i, &dU_center, &fl, &fr]
	(int kl, int ku, int jl, int ju, int il, int iu) {
	for (int iz=kl; iz<ku; iz++) {
	  for (int iy=jl; iy<ju; iy++) {
	    for (int ix=il; ix<iu; ix++) {
	      dU_center(iz,iy,ix) -= dtdx_i * (fr(iz,iy,ix) - fl(iz,iy,ix));
	    }
	  }
	}
      };

      if (masked){
	mask_iter(mask, inner_loop, dim, 3);
      } else {
	inner_loop(0, dU_center.shape(0), 0, dU_center.shape(1),
		   0, dU_center.shape(2));
      }
    }
  }
}

//----------------------------------------------------------------------

bool EnzoIntegrableUpdate::identify_fallback
(Block *block, Grouping &initial_integrable_group, Grouping &dUcons_group,
 CelloArray<bool,3> &fallback_mask, int stale_depth) const
{
  // In the future, this should prepare a mask indicating exactly which
  // locations require recomputed values of dU relying on the fallback solver
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EFlt3DArray rho, drho;
  rho = array_factory.from_grouping(initial_integrable_group, "density", 0);
  drho = array_factory.from_grouping(dUcons_group, "density", 0);

  CSlice slc(stale_depth, -1*stale_depth);
  CelloArray<bool,3> mask = fallback_mask.subarray(slc, slc, slc);
  ASSERT("EnzoIntegrableUpdate::identify_fallback",
	 "The density field and fallback_mask must have the same shape",
	 (rho.shape(0) == mask.shape(0) && rho.shape(1) == mask.shape(1) &&
	  rho.shape(2) == mask.shape(2)));
  bool out = false;
  for (int iz=0; iz<rho.shape(0); iz++) {
    for (int iy=0; iy<rho.shape(1); iy++) {
      for (int ix=0; ix<rho.shape(2); ix++) {
	mask(iz,iy,ix) = rho(iz,iy,ix) < (-1*drho(iz,iy,ix));
	if ( mask(iz,iy,ix) ) {
	  WARNING5("EnzoIntegrableUpdate::identify_fallback",
		   ("At (iz,iy,ix) = (%d,%d,%d), the total change in "
		    "density, %.15e, causes density, %.15e, to go negative\n"),
		   iz,iy,ix,drho(iz,iy,ix),rho(iz,iy,ix));
	  out = true;
	}
      }
    }
  }
  return out;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_quantities
(Block *block, Grouping &initial_integrable_group, Grouping &dUcons_group,
 Grouping &out_integrable_group, Grouping &out_conserved_passive_scalar,
 EnzoEquationOfState *eos, int stale_depth) const
{

  // Update passive scalars, it doesn't currently support renormalizing to 1
  update_passive_scalars_(block, initial_integrable_group, dUcons_group,
			  out_conserved_passive_scalar, stale_depth);

  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos->get_density_floor();

  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *cur_prim, *dU, *out_prim;
  cur_prim = registry.load_array_of_fields(block, lut_, nfields_,
					   initial_integrable_group,
					   stale_depth);
  dU       = registry.load_array_of_fields(block, lut_, nfields_,
					   dUcons_group, stale_depth);
  out_prim = registry.load_array_of_fields(block, lut_, nfields_,
					   out_integrable_group, stale_depth);

  for (int iz = 1; iz < (cur_prim[lut_.density].shape(0) - 1); iz++) {
    for (int iy = 1; iy < (cur_prim[lut_.density].shape(1) - 1); iy++) {
      for (int ix = 1; ix < (cur_prim[lut_.density].shape(2) - 1); ix++) {

	// get the initial density
	enzo_float old_rho = cur_prim[lut_.density](iz,iy,ix);

	// now update the integrable primitives that are conserved
	for (int i = conserved_start_; i < conserved_stop_; i++){
	  out_prim[i](iz,iy,ix) = cur_prim[i](iz,iy,ix) + dU[i](iz,iy,ix);
	}

	// possibly place a floor on the updated density.
	enzo_float new_rho = out_prim[lut_.density](iz,iy,ix);
	new_rho = EnzoEquationOfState::apply_floor(new_rho, density_floor);
	out_prim[lut_.density](iz,iy,ix) = new_rho;

	// update the specific integrable primitives
	enzo_float inv_new_rho = 1./new_rho;
	for (int i = specific_start_; i < specific_stop_; i++){
	  out_prim[i](iz,iy,ix) =
	    (cur_prim[i](iz,iy,ix) * old_rho + dU[i](iz,iy,ix)) * inv_new_rho;
	}

      }
    }
  }

  // apply floor to energy and sync the internal energy with total energy
  // (the latter only occurs if the dual energy formalism is in use)
  eos->apply_floor_to_energy_and_sync(block, out_integrable_group,
				      stale_depth+1);

  delete[] cur_prim;  delete[] dU;  delete[] out_prim;
}

//----------------------------------------------------------------------

void EnzoIntegrableUpdate::update_passive_scalars_
  (Block *block, Grouping &initial_integrable_group, Grouping &dUcons_group,
   Grouping &out_conserved_passive_scalar, int stale_depth) const
{

  // Does not currently handle passively advected scalars that must sum to 1
  // The easiest way to do that would be to divide each group of fields that
  // must sum to 1 into different groupings

  EnzoFieldArrayFactory array_factory(block, stale_depth);

  EFlt3DArray cur_rho
    = array_factory.from_grouping(initial_integrable_group, "density", 0);

  // cell-centered grid dimensions
  int mz, my, mx;
  mz = cur_rho.shape(0);    my = cur_rho.shape(1);    mx = cur_rho.shape(2);

  std::vector<std::string> group_names = this->passive_groups_;
  for (std::string group_name : group_names){
    int num_fields = initial_integrable_group.size(group_name);

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      EFlt3DArray cur_specific, out_conserved, dU;
      cur_specific = array_factory.from_grouping(initial_integrable_group,
						 group_name, field_ind);
      out_conserved = array_factory.from_grouping(out_conserved_passive_scalar,
						  group_name, field_ind);
      dU = array_factory.from_grouping(dUcons_group, group_name, field_ind);

      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {

	    out_conserved(iz,iy,ix)
	      = (cur_specific(iz,iy,ix) * cur_rho(iz,iy,ix) + dU(iz,iy,ix));
	  }
	}
      }

    }
  }
}
