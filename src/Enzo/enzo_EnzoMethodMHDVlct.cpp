// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDVlct.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri June 14 2019
/// @brief    [\ref Enzo] Implementation of the EnzoMethodMHDVlct class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"
#include <algorithm>    // std::copy

//----------------------------------------------------------------------

EnzoMethodMHDVlct::EnzoMethodMHDVlct (std::string rsolver,
				      std::string half_recon_name,
				      std::string full_recon_name,
				      double gamma, double theta_limiter,
				      double density_floor,
				      double pressure_floor,
				      std::string mhd_choice,
				      bool dual_energy_formalism,
				      double dual_energy_formalism_eta)
  : Method()
{
  // Initialize equation of state (check the validity of quantity floors)
  EnzoEquationOfState::check_floor(density_floor);
  EnzoEquationOfState::check_floor(pressure_floor);
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor,
			  dual_energy_formalism, dual_energy_formalism_eta);

  // Determine whether magnetic fields are to be used
  mhd_choice_ = parse_bfield_choice_(mhd_choice);

  // determine integrable and reconstructable quantities (and passive scalars)
  determine_quantities_(eos_, integrable_group_names_,
			reconstructable_group_names_, passive_group_names_);

  FieldDescr * field_descr = cello::field_descr();

  // "pressure" is only used to compute the timestep
  ASSERT("EnzoMethodMHDVlct", "\"pressure\" must be a permanent field",
	 field_descr->is_field("pressure"));


  // setup primitive_group_, and bfieldi_group_
  // (also checks that the integrable fields of primitive_group_ and all the
  // fields of bfieldi_group_ exist and are permanent)
  setup_groupings_(integrable_group_names_, reconstructable_group_names_,
		   passive_group_names_);


  // Initialize the default Refresh object - May want to adjust
  // number of ghost zones based on reconstructor choice.
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // Add cell-centered fields holding actively advected integrable quantities
  add_group_fields_to_refresh_(refresh, *primitive_group_,
                               integrable_group_names_);
  // Add cell-centered fields holding the conserved form of the passively
  // advected quantities
  add_group_fields_to_refresh_(refresh, *(field_descr->groups()),
			       passive_group_names_);

  /// Add interface fields (if necessary) to the refresh list
  if (mhd_choice_ == bfield_choice::constrained_transport) {
    EnzoConstrainedTransport::update_refresh(refresh);
  }

  // Initialize the remaining component objects
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, half_recon_name, (enzo_float)theta_limiter);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, full_recon_name, (enzo_float)theta_limiter);
  riemann_solver_ = EnzoRiemann::construct_riemann
    (integrable_group_names_,      passive_group_names_, rsolver);
  integrable_updater_ = new EnzoIntegrableUpdate(integrable_group_names_,
						 true, passive_group_names_);
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::bfield_choice EnzoMethodMHDVlct::parse_bfield_choice_
(std::string choice) const noexcept
{
  std::string formatted(choice.size(), ' ');
  std::transform(choice.begin(), choice.end(), formatted.begin(),
		 ::tolower);
  if (formatted == std::string("no_bfield")){
    return bfield_choice::no_bfield;
  } else if (formatted == std::string("unsafe_constant_uniform")){
    ERROR("EnzoMethodMHDVlct::parse_bfield_choice_",
          "constant_uniform is primarilly for debugging purposes. DON'T use "
          "for science runs (things can break).");
    return bfield_choice::unsafe_const_uniform;
  } else if (formatted == std::string("constrained_transport")){
    return bfield_choice::constrained_transport;
  } else {
    ERROR("EnzoMethodMHDVlct::parse_bfield_choice_",
          "Unrecognized choice. Known options include \"no_bfield\" and "
          "\"constrained_transport\"");
    return bfield_choice::no_bfield;
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::determine_quantities_
(EnzoEquationOfState *eos, std::vector<std::string> &integrable_quantities,
 std::vector<std::string> &reconstructable_quantities,
 std::vector<std::string> &passive_groups)
{
#ifdef CONFIG_USE_GRACKLE
  if (enzo::config()->method_grackle_use_grackle){
    // make sure all the required fields are defined so that the group of
    // "colour" fields is accurate (needed for identifying passive scalars)
    EnzoMethodGrackle::define_required_grackle_fields();
    // Not quite ready to support a variable gamma
  }
#endif


  std::string common[] {"density", "velocity"};
  for (std::string quantity : common){
    integrable_quantities.push_back(quantity);
    reconstructable_quantities.push_back(quantity);
  }

  if (mhd_choice_ != bfield_choice::no_bfield){
    integrable_quantities.push_back("bfield");
    reconstructable_quantities.push_back("bfield");
  }

  if (eos->is_barotropic()){
    ERROR("EnzoMethodMHDVlct::determine_quantities_",
	  "Not presently equipped to handle barotropic equations of state.");
  } else {
    // add specific total energy to integrable quantities
    integrable_quantities.push_back("total_energy");
    // add pressure to reconstructable quantities
    reconstructable_quantities.push_back("pressure");
    // add specific internal energy to integrable quantities (if using the dual
    // engery formalism)
    if (eos->uses_dual_energy_formalism()){
      integrable_quantities.push_back("internal_energy");
    }
  }

  // Now to setup the list of passively advected group names
  // retrieve list of all possible group names of passively advected scalars
  std::vector<std::string> all_passive_group_names =
    EnzoCenteredFieldRegistry::passive_scalar_group_names();

  // We'll now iterate through all known possible groups names that could
  // contain passive scalars and check to see if any fields were labelled as
  // being a part of these groups in the input file (if so, then add the group
  // name to the list of passively advected group names)
  Grouping* reference_grouping = cello::field_descr()->groups();
  for (std::size_t i=0; i < all_passive_group_names.size(); i++){
    std::string group_name = all_passive_group_names[i];

    if (reference_grouping->size(group_name) > 0){
      passive_groups.push_back(group_name);
    }
  }

}

//----------------------------------------------------------------------

void add_passive_groups_(Grouping &grouping, const std::string field_prefix,
			 const std::vector<std::string> group_names)
{
  // Helper function that adds entries of passively advected groups (specified
  // in the configuration file) to an existing Grouping object

  // This includes the groups specified in the parameter file
  Grouping* reference_grouping = cello::field_descr()->groups();

  for (unsigned int i=0;i<group_names.size();i++){

    std::string group_name = group_names[i];
    int num_fields = reference_grouping->size(group_name);

    for (int j=0; j < num_fields; j++){
      // Determine field_name
      std::string field_name = (field_prefix +
				reference_grouping->item(group_name,j));

      // add the field to the grouping
      grouping.add(field_name, group_name);
    }
  }
}

//----------------------------------------------------------------------

// Returns the unique members of a combination of 2 vectors
std::vector<std::string> unique_combination_(const std::vector<std::string> &a,
					     const std::vector<std::string> &b)
{
  std::vector<std::string> out;
  // copy elements from a
  for (std::size_t i = 0; i<a.size(); i++){
    std::string name = a[i];
    if (std::find(out.begin(), out.end(), name) == out.end()){
      out.push_back(name);
    }
  }
  // copy elements from b
  for (std::size_t i = 0; i<b.size(); i++){
    std::string name = b[i];
    if (std::find(out.begin(), out.end(), name) == out.end()){
      out.push_back(name);
    }
  }
  return out;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::setup_groupings_
(std::vector<std::string> &integrable_groups,
 std::vector<std::string> &reconstructable_groups,
 std::vector<std::string> &passive_groups)
{

  FieldDescr * field_descr = cello::field_descr();

  // first come up with a vector group names that represents the union of
  // integrable_groups and reconstructable_groups
  std::vector<std::string> groups;
  groups = unique_combination_(integrable_groups,reconstructable_groups);

  // now setup primitive_group_ using the names in groups
  primitive_group_ = EnzoCenteredFieldRegistry::build_grouping(groups, "");

  // We should check that all the fields in integrable groups are real
  // permenant fields
  for (std::size_t i = 0; i<integrable_groups.size(); i++){
    std::string group_name = integrable_groups[i];
    int num_fields = primitive_group_->size(group_name);

    for (int j = 0; j<num_fields; j++){
      std::string field_name = primitive_group_->item(group_name,j);

      ASSERT1("EnzoMethodMHDVlct::setup_groupings_",
	      "\"%s\" must be the name of a permanent field",
	      field_name.c_str(), field_descr->is_field(field_name));
    }
  }

  // Add temporary fields to primitive_group_ to hold specific forms (mass
  // fraction) of passively advected scalars. These get are used to hold the
  // values after they are converted from conserved form (density) which are
  // tracked in permanent fields
  add_passive_groups_(*primitive_group_, "specific_passive_", passive_groups);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::add_group_fields_to_refresh_
(Refresh * refresh, Grouping &grouping, std::vector<std::string> group_names)
{
  FieldDescr * field_descr = cello::field_descr();

  for (unsigned int i=0;i<group_names.size();i++){

    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);

    for (int j=0; j < num_fields; j++){
      // Determine field_name
      std::string field_name = grouping.item(group_name,j);
      ASSERT1("EnzoMethodMHDVlct::add_group_fields_to_refresh_",
	      "%s must be a permanent field",field_name.c_str(),
	      field_descr->is_permanent(field_name))

      // add the field to the refresh object
      refresh->add_field(field_descr->field_id(field_name));
    }
  }
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete primitive_group_;

  delete eos_;
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
  const bool up = p.isUnpacking();

  p|eos_;
  p|integrable_group_names_;
  p|reconstructable_group_names_;
  p|passive_group_names_;

  int has_prim_group = (primitive_group_ != nullptr);
  p|has_prim_group;
  if (has_prim_group){
    if (up){
      primitive_group_ = new Grouping;
    }
    p|*primitive_group_;
  } else {
    primitive_group_ = nullptr;
  }

  // sanity check:
  ASSERT("EnzoMethodMHDVlct::pup", "primitive_group_ should not be NULL",
         primitive_group_ != nullptr);

  p|half_dt_recon_;
  p|full_dt_recon_;
  p|riemann_solver_;
  p|integrable_updater_;
  p|mhd_choice_;
}

//----------------------------------------------------------------------

void add_arrays_to_map_(Block * block,
                        Grouping& grouping,
                        const std::vector<std::string>& group_names,
                        int dim, EnzoEFltArrayMap& map,
                        bool allow_missing_group, bool enforce_num_Groups,
                        Grouping* ref_grouping)
{

  char suffixes[3] = {'x','y','z'};
  EnzoFieldArrayFactory array_factory(block,0);

  for (std::string group_name : group_names){
    int num_fields = grouping.size(group_name);

    if ((num_fields == 0) && allow_missing_group){ continue; }

    ASSERT("EnzoMethodMHDVlct::compute_flux_",
           "all groups must have 1 or 3 fields.",
           !enforce_num_Groups || ((num_fields == 1) || (num_fields == 3)));

    for (int field_ind=0; field_ind<num_fields; field_ind++){
      std::string field_name = grouping.item(group_name,field_ind);

      std::string key;
      if (ref_grouping != nullptr){
        key = ref_grouping->item(group_name, field_ind);
      } else if (num_fields == 3){
        key = group_name;
        key.push_back('_');
        key.push_back(suffixes[field_ind]);
      } else {
        key = group_name;
      }

      if (map.contains(key)){
        ERROR1("EnzoEFltArrayMap::from_grouping",
               "EnzoEFltArrayMap can't hold more than one field called \"%s\"",
               key.c_str());
      }

      if (dim == -1){
        map[key] = array_factory.from_name(field_name);
      } else {
        map[key] = array_factory.assigned_center_from_name(field_name, dim);
      }

    }

  }
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::nonpassive_primitive_map_(Block * block)
  const throw ()
{
  EnzoEFltArrayMap primitive_map("primitive");
  // need to combine group_names to handle the internal energy source term
  std::vector<std::string> all_prim_group_names =
    unique_combination_(reconstructable_group_names_,
                        integrable_group_names_);
  add_arrays_to_map_(block, *primitive_group_, all_prim_group_names,
                     -1, primitive_map, false, true, NULL);
  return primitive_map;
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::conserved_passive_scalar_map_
(Block * block) const throw ()
{
  // get Grouping of field names that store passively advected scalars in
  // conserved-form
  Grouping *conserved_passive_scalars = cello::field_descr()->groups();
  EnzoEFltArrayMap conserved_passive_scalar_map("conserved_passive_scalar");
  add_arrays_to_map_(block, *conserved_passive_scalars,
                     passive_group_names_, -1,
                     conserved_passive_scalar_map, false, false,
                     conserved_passive_scalars);
  return conserved_passive_scalar_map;
}

void update_flux_data_(Block * block, const EnzoEFltArrayMap &flux_map,
                       int dim)
{
  Field field = block->data()->field();
  const EnzoPermutedCoordinates coord(dim);

  // first get dimensions of cell-centered fields along i,j,k
  int ghost[3], cc_shape[3]; // the values are ordered as x,y,z
  field.ghost_depth(0, &(ghost[0]), &(ghost[1]), &(ghost[2]));
  field.dimensions(0, &(cc_shape[0]), &(cc_shape[1]), &(cc_shape[2]));

  int gi = ghost[coord.i_axis()];     int mi = cc_shape[coord.i_axis()];
  int gj = ghost[coord.j_axis()];     int mj = cc_shape[coord.j_axis()];
  int gk = ghost[coord.k_axis()];     int mk = cc_shape[coord.k_axis()];

  // now, compute slices along j and k axes dimension. We build these with
  // exact indices (i.e. we don't use negative indices) just in case the arrays
  // holding the fluxes are bigger than they need to be
  const CSlice j_slc(gj, mj - gj); // slices along j and k just include the
  const CSlice k_slc(gk, mk - gk); // active zone

  // the fluxes on the lower face of the block are found between the cells at
  // (gi-1) and gi. Because the flux arrays have space for values located on
  // all interior faces between cells, these values are located at the index
  // (gi-1).
  const CSlice i_lower_slc(gi - 1, gi);
  // the fluxes on the upper face of the block are found between the cells at
  // (mi-gi-1) and (mi-gi). Because the flux arrays have space for values
  // located on all interior faces between cells, these values are located at
  // the index (mi-gi-1).
  const CSlice i_upper_slc(mi - gi - 1, mi -gi);

  // Finally we actually store the fluxes
  FluxData * flux_data = block->data()->flux_data();

  const int nf = flux_data->num_fields();
  for (int i_f=0; i_f <nf; i_f++) {
    int * flux_index = 0;
    const int index_field = flux_data->index_field(i_f);
    const std::string field_name = field.field_name(index_field);

    // note the field_name is the same as the key

    for (int face = 0; face < 2; face++){

      FaceFluxes * ff_b = flux_data->block_fluxes(dim,face,i_f);
      int mx, my, mz;
      ff_b->get_size(&mx,&my,&mz);
      int dx, dy, dz;
      enzo_float* dest = ff_b->flux_array(&dx,&dy,&dz).data();

      const CSlice& i_slc = (face == 0) ? i_lower_slc : i_upper_slc;
      EFlt3DArray tmp = flux_map.at(field_name);
      EFlt3DArray src = coord.get_subarray(tmp, k_slc, j_slc, i_slc);

      // as a sanity check: may want to check
      ASSERT7("update_flux_data_",
              ("For the flux along %d, the subarray is expected to have shape "
               "(%d,%d,%d). Instead, its shape is (%d,%d,%d)."),
              dim, mz,my,mx, src.shape(0), src.shape(1), src.shape(2),
              ((src.shape(0) == mz) && (src.shape(1) == my) &&
               (src.shape(2) == mx))
              );

      for (int iz = 0; iz < mz; iz++){
        for (int iy = 0; iy < my; iy++){
          for (int ix = 0; ix < mx; ix++){
            dest[iz*dz + iy*dy + ix*dx] = src(iz,iy,ix);
          }
        }
      }

    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  // the following is getting copied from EnzoMethodPpm::compute
  Field field = block->data()->field();
  auto field_names = field.groups()->group_list("conserved");
  const int nf = field_names.size();
  std::vector<int> field_list;
  field_list.resize(nf);
  for (int i=0; i<nf; i++) {
    field_list[i] = field.field_id(field_names[i]);
  }

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  // do we need to allocate every cycle?
  block->data()->flux_data()->allocate (nx,ny,nz,field_list);

  if (block->is_leaf()) {
    // Check that the mesh size and ghost depths are appropriate
    check_mesh_and_ghost_size_(block);

    // declaring Maps of arrays and stand-alone arrays that wrap existing
    // fields and/or serve as scratch space.

    // map that holds arrays wrapping the Cello Fields holding each of the
    // primitive quantities. Additionally, this also includes temporary arrays
    // used to hold the specific form of the passive scalar
    EnzoEFltArrayMap primitive_map; // this will be overwritten

    // map used for storing primitive values at the half time-step. This
    // includes key,array pairs for each entry in primitive_map (there should
    // be no aliased fields shared between maps)
    EnzoEFltArrayMap temp_primitive_map("temp_primitive");

    // map holding the arrays wrapping each fields corresponding to a passively
    // advected scalar. The scalar is in conserved form.
    EnzoEFltArrayMap conserved_passive_scalar_map; // this will be overwritten

    // holds left and right reconstructed primitives (scratch-space)
    EnzoEFltArrayMap priml_map("priml");
    EnzoEFltArrayMap primr_map("primr");

    // Arrays used to store the pressure computed from the reconstructed left
    // and right primitives
    // Note: in the case of adiabatic fluids, pressure is a reconstructable
    //       quantity and entries are included for it in priml_map and
    //       primr_map. In that case, pressure_l and pressure_r are aliases of
    //       those arrays.
    EFlt3DArray pressure_l, pressure_r;

    // maps used to store fluxes (in the future, these will wrap FluxData
    // entries)
    EnzoEFltArrayMap xflux_map("xflux");
    EnzoEFltArrayMap yflux_map("yflux");
    EnzoEFltArrayMap zflux_map("zflux");

    // map of arrays  used to accumulate the changes to the conserved forms of
    // the integrable quantities and passively advected scalars. In other
    // words, at the start of the (partial) timestep, the fields are all set to
    // zero and are used to accumulate the flux divergence and source terms. If
    // CT is used, it won't have space to store changes in the magnetic fields.
    EnzoEFltArrayMap dUcons_map("dUcons");

    // This is a list of lists of passive scalar names (or keys). The first
    // sublist holds the names of all quantities that undergo normal passive
    // advection. Subsequent lists hold sets of names for scalars whose values
    // must sum to 1.0 (like species).
    std::vector<std::vector<std::string>> passive_lists {{}};

    setup_arrays_(block, primitive_map, temp_primitive_map,
                  conserved_passive_scalar_map, priml_map, primr_map,
		  pressure_l, pressure_r, xflux_map, yflux_map, zflux_map,
                  dUcons_map, passive_lists);

    // Setup a pointer to an array that used to store interface velocity fields
    // from computed by the Riemann Solver (to use in the calculation of the
    // internal energy source term). If the dual energy formalism is not in
    // use, don't actually allocate the array and set the pointer to NULL.
    EFlt3DArray interface_velocity_arr, *interface_velocity_arr_ptr;
    if (eos_->uses_dual_energy_formalism()){
      EFlt3DArray density = primitive_map.at("density");
      interface_velocity_arr = EFlt3DArray(density.shape(0), density.shape(1),
                                           density.shape(2));
      interface_velocity_arr_ptr = &interface_velocity_arr;
    } else {
      interface_velocity_arr_ptr = NULL;
    }

    // allocate constrained transport object
    EnzoConstrainedTransport *ct = NULL;
    if (mhd_choice_ == bfield_choice::constrained_transport) {
      ct = new EnzoConstrainedTransport(block, 2);
    }

    const double* const cell_widths = enzo::block(block)->CellWidth;

    double dt = block->dt();

    // stale_depth indicates the number of field entries from the outermost
    // field value that the region including "stale" values (need to be
    // refreshed) extends over.
    int stale_depth = 0;

    // convert the passive scalars from conserved form to specific form
    // (outside the integrator, they are treated like conserved densities)
    compute_specific_passive_scalars_(passive_lists, primitive_map["density"],
                                      conserved_passive_scalar_map,
                                      primitive_map, stale_depth);

    // repeat the following loop twice (for half time-step and full time-step)

    for (int i=0;i<2;i++){
      double cur_dt = (i == 0) ? dt/2. : dt;
      EnzoEFltArrayMap& cur_integrable_map =
        (i == 0) ? primitive_map      : temp_primitive_map;
      EnzoEFltArrayMap& out_integrable_map =
        (i == 0) ? temp_primitive_map :      primitive_map;
      // For the purposes of making the calculation procedure slightly more
      // more explicit, we distinguish between cur_integrable_group and
      // cur_reconstructable_group. Due to the high level of overlap between
      // these, they are simply aliases of the same underlying Grouping that
      // holds groups for both of them
      EnzoEFltArrayMap& cur_reconstructable_map = cur_integrable_map;

      EnzoReconstructor *reconstructor;

      if (i == 0){
        reconstructor = half_dt_recon_;
        // ct does NOT need to be incremented
      } else {
        reconstructor = full_dt_recon_;

        if (ct != NULL){ ct->increment_partial_timestep(); }

        // After the fluxes were added to the passive scalar in the first half
        // timestep, the values were stored in conserved form in the fields
        // held by conserved_passive_scalar_map.
        // Need to convert them to specific form
        compute_specific_passive_scalars_(passive_lists,
                                          cur_integrable_map["density"],
                                          conserved_passive_scalar_map,
                                          cur_integrable_map, stale_depth);
      }

      // set all elements of the arrays in dUcons_map to 0 (throughout the rest
      // of the current loop, flux divergence and source terms will be
      // accumulated in these arrays)
      integrable_updater_->clear_dUcons_map(dUcons_map, 0., passive_lists);

      // Compute the reconstructable quantities from the integrable quantites
      // Although cur_integrable_map holds the passive scalars in integrable
      // form, the conserved form of the values is required in case Grackle is
      // being used.
      //
      // Note: cur_integrable_map and cur_reconstructable_map are aliases of
      // the same map since there is such a large degree of overlap between
      // reconstructable and integrable quantities
      //
      // For a barotropic gas, the following nominally does nothing
      // For a non-barotropic gas, the following nominally computes pressure
      eos_->reconstructable_from_integrable(cur_integrable_map,
                                            cur_reconstructable_map,
                                            conserved_passive_scalar_map,
                                            stale_depth, passive_lists);

      // Compute flux along each dimension
      compute_flux_(0, cur_dt, cell_widths[0], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    xflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, ct, stale_depth, passive_lists);
      compute_flux_(1, cur_dt, cell_widths[1], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    yflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, ct, stale_depth, passive_lists);
      compute_flux_(2, cur_dt, cell_widths[2], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    zflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, ct, stale_depth, passive_lists);

      if (i == 1 && ct == NULL) {
        if (eos_->uses_dual_energy_formalism()){
          // the interface velocity on the edge of the block will be different!
          // I think the answer here is to add this effect to the internal
          // energy flux
          ERROR("EnzoMethodMHDVlct::compute",
                "Handling of the dual energy source term is won't be properly "
                "corrected!");
        }
        update_flux_data_(block, xflux_map, 0);
        update_flux_data_(block, yflux_map, 1);
        update_flux_data_(block, zflux_map, 2);
      }

      // increment the stale_depth
      stale_depth+=reconstructor->immediate_staling_rate();

      // This is where source terms should be computed (added to dUcons_group)

      // Update Bfields
      if (ct != NULL) {
        ct->update_all_bfield_components(cur_integrable_map, xflux_map,
                                         yflux_map, zflux_map,
                                         out_integrable_map, cur_dt,
                                         stale_depth);
      }

      // Update quantities (includes flux divergence and source terms) 
      // This needs to happen after updating the cell-centered B-field so that
      // the pressure floor can be applied to the total energy (and if
      // necessary the total energy can be synchronized with internal energy)
      //
      // Note: updated passive scalars are NOT saved in out_integrable_group in
      //     specific form. Instead they are saved inconserved_passive_scalars
      //     in conserved form.
      integrable_updater_->update_quantities(primitive_map, dUcons_map,
                                             out_integrable_map,
                                             conserved_passive_scalar_map,
                                             eos_, stale_depth, passive_lists);

      // increment stale_depth since the inner values have been updated
      // but the outer values have not
      stale_depth+=reconstructor->delayed_staling_rate();
    }

    if (ct != NULL){delete ct;}
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::check_mesh_and_ghost_size_(Block *block) const
{
  Field field = block->data()->field();

  // Check that the mesh size is appropriate
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
  ASSERT("EnzoMethodMHDVlct::compute",
	 "Active zones on each block must be >= ghost depth.",
	 nx >= gx && ny >= gy && nz >= gz);

  // Check that the (cell-centered) ghost depth is large enough
  // Face-centered ghost could in principle be 1 smaller
  int min_ghost_depth = (half_dt_recon_->total_staling_rate() +
			 full_dt_recon_->total_staling_rate());
  ASSERT1("EnzoMethodMHDVlct::compute", "ghost depth must be at least %d.",
	  min_ghost_depth, std::min(nx, std::min(ny, nz)) >= min_ghost_depth);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_specific_passive_scalars_
(const std::vector<std::vector<std::string>> passive_lists,
 EFlt3DArray& density, EnzoEFltArrayMap& conserved_passive_scalar_map,
 EnzoEFltArrayMap& specific_passive_scalar_map, int stale_depth) const noexcept
{
  int mz = density.shape(0);
  int my = density.shape(1);
  int mx = density.shape(2);

  for (const std::vector<std::string> cur_l : passive_lists){
    for (const std::string& key : cur_l){

      EFlt3DArray cur_conserved = conserved_passive_scalar_map.at(key);
      EFlt3DArray out_specific = specific_passive_scalar_map.at(key);

      for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
	for (int iy = stale_depth; iy < my - stale_depth; iy++) {
	  for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
	    out_specific(iz,iy,ix) = cur_conserved(iz,iy,ix)/density(iz,iy,ix);
	  }
	}
      }

    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_flux_
(int dim, double cur_dt, enzo_float cell_width,
 EnzoEFltArrayMap &reconstructable_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
 EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
 EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
 EnzoConstrainedTransport *ct_handler, int stale_depth,
 const std::vector<std::vector<std::string>>& passive_lists) const noexcept
{

  // purely for the purposes of making the caluclation more explicit, we define
  // the following aliases for priml_map and primr_map
  EnzoEFltArrayMap &reconstructable_l = priml_map;
  EnzoEFltArrayMap &reconstructable_r = primr_map;
  EnzoEFltArrayMap &integrable_l = priml_map;
  EnzoEFltArrayMap &integrable_r = primr_map;

  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(reconstructable_map,
                                      reconstructable_l, reconstructable_r,
				      dim, eos_, stale_depth,
                                      passive_lists);

  // We temporarily increment the stale_depth for the rest of this calculation
  // here. We can't fully increment otherwise it will screw up the
  // reconstruction along other dimensions
  int cur_stale_depth = stale_depth + reconstructor.immediate_staling_rate();

  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces (should potentially be handled internally by reconstructor)
  if (ct_handler != NULL) {
    ct_handler->correct_reconstructed_bfield(reconstructable_l,
                                             reconstructable_r,
					     dim, cur_stale_depth);
  }

  // Calculate integrable values on left and right faces:
  eos_->integrable_from_reconstructable(reconstructable_l, integrable_l,
					cur_stale_depth, passive_lists);
  eos_->integrable_from_reconstructable(reconstructable_r, integrable_r,
					cur_stale_depth, passive_lists);

  // Calculate pressure on left and right faces:
  eos_->pressure_from_reconstructable(reconstructable_l, pressure_l,
                                      cur_stale_depth);
  eos_->pressure_from_reconstructable(reconstructable_r, pressure_r,
                                      cur_stale_depth);

  // Next, compute the fluxes
  riemann_solver_->solve(integrable_l, integrable_r, pressure_l, pressure_r,
                         flux_map, dim, eos_, cur_stale_depth, passive_lists,
                         interface_velocity_arr_ptr);

  // Accumulate the change in integrable quantities from these flux_map in
  // dUcons_map
  integrable_updater_->accumulate_flux_component(dim, cur_dt, cell_width,
                                                 flux_map, dUcons_map,
                                                 cur_stale_depth,
                                                 passive_lists);

  // if using dual energy formalism, compute the dual energy formalism
  if (eos_->uses_dual_energy_formalism()){
    EnzoSourceInternalEnergy eint_src;
    eint_src.calculate_source(dim, cur_dt, cell_width, reconstructable_map,
                              dUcons_map, *interface_velocity_arr_ptr, eos_,
                              cur_stale_depth);
  }

  // Finally, have the ct handler record the upwind direction
  if (ct_handler != NULL){
    ct_handler->identify_upwind(flux_map, dim, cur_stale_depth);
  }
}

//----------------------------------------------------------------------

void add_temporary_arrays_to_map_
(EnzoEFltArrayMap &map, std::array<int,3> &shape,
 const std::vector<std::string>* const names,
 const std::vector<std::vector<std::string>>* const passive_lists,
 bool skip_unregisterred_names = false)
{

  if (names != nullptr){
    for (const std::string& name : (*names)){
      bool success, is_vector;
      success = EnzoCenteredFieldRegistry::quantity_properties (name,
                                                                &is_vector);
      if (skip_unregisterred_names & (!success)){
        continue;
      } else {
        ASSERT1("add_temporary_arrays_to_map_",
                ("\"%s\" is not registered in EnzoCenteredFieldRegistry"),
                name.c_str(), success);
      }

      if (is_vector){
        map[name + "_x"] = EFlt3DArray(shape[0], shape[1], shape[2]);
        map[name + "_y"] = EFlt3DArray(shape[0], shape[1], shape[2]);
        map[name + "_z"] = EFlt3DArray(shape[0], shape[1], shape[2]);
      } else {
        map[name] = EFlt3DArray(shape[0], shape[1], shape[2]);
      }
    }
  }

  if (passive_lists != nullptr){
    for (const std::vector<std::string>& cur_list : (*passive_lists)){
      for (const std::string& key : cur_list){
        map[key] = EFlt3DArray(shape[0], shape[1], shape[2]);
      }
    }
  }

}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::setup_arrays_
(Block *block, EnzoEFltArrayMap &primitive_map,
 EnzoEFltArrayMap &temp_primitive_map,
 EnzoEFltArrayMap &conserved_passive_scalar_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
 EnzoEFltArrayMap &xflux_map, EnzoEFltArrayMap &yflux_map,
 EnzoEFltArrayMap &zflux_map, EnzoEFltArrayMap &dUcons_map,
 std::vector<std::vector<std::string>> &passive_lists) noexcept
{

  // allocate stuff! Make sure to do it in a way such that we don't have to
  // separately deallocate it!

  // setup passive_lists. It's a list of lists of passive scalar names. The
  // first sublist holds all names that are normally passively advected.
  // Subsequent lists hold sets of names for scalars whose values must sum to
  // 1 (like species).
  Grouping *conserved_passive_scalar_grouping = cello::field_descr()->groups();
  for (std::string group_name : passive_group_names_){
    int num_fields = conserved_passive_scalar_grouping->size(group_name);
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      std::string field_name =
        conserved_passive_scalar_grouping->item(group_name, field_ind);
      passive_lists[0].push_back(field_name);
    }
  }

  // To assist with setting up arrays, let's create a list of ALL primitive
  // keys (including passive scalars) and all integrable keys. These are the
  // same thing except the latter excludes quantities (like pressure) that we
  // don't compute pressure for.
  std::vector<std::string> combined_key_list = unique_combination_
    (reconstructable_group_names_, integrable_group_names_);

  // First, setup conserved_passive_scalar_map
  conserved_passive_scalar_map = conserved_passive_scalar_map_(block);

  // Next, setup nonpassive components of primitive_map
  primitive_map = nonpassive_primitive_map_(block);
  std::array<int,3> shape = {primitive_map.at("density").shape(0),
                             primitive_map.at("density").shape(1),
                             primitive_map.at("density").shape(2)};

  add_temporary_arrays_to_map_(primitive_map, shape, nullptr, &passive_lists);

  // Then, setup temp_primitive_map
  add_temporary_arrays_to_map_(temp_primitive_map, shape, &combined_key_list,
                               &passive_lists);

  // Prepare arrays to hold fluxes (it should include groups for all actively
  // and passively advected quantities)
  EnzoEFltArrayMap* flux_maps[3] = {&zflux_map, &yflux_map, &xflux_map};
  for (std::size_t i = 0; i < 3; i++){
    std::array<int,3> cur_shape = shape; // makes a deep copy
    cur_shape[i] -= 1;
    add_temporary_arrays_to_map_(*(flux_maps[i]), cur_shape,
                                 &combined_key_list, &passive_lists);
  }

  // Prepare fields used to accumulate all changes to the actively advected and
  // passively advected quantities. If CT is in use, dUcons_group should not
  // have storage for magnetic fields since CT independently updates magnetic
  // fields (this exclusion is implicitly handled integrable_updater_)
  const std::vector<std::string> temp_l =
    integrable_updater_->combined_integrable_groups();
  add_temporary_arrays_to_map_(dUcons_map, shape, &temp_l, &passive_lists,
                               true);

  // Prepare temporary fields for priml and primr
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  add_temporary_arrays_to_map_(priml_map, shape, &combined_key_list,
                               &passive_lists);
  add_temporary_arrays_to_map_(primr_map, shape, &combined_key_list,
                               &passive_lists);

  // If there are pressure entries in priml_map and primr_map (depends on the
  // EOS), set pressure_l and pressure_name_r equal to
  // those field names. Otherwise, reserve/allocate left/right pressure fields
  if (priml_map.contains("pressure")) {
    pressure_l = priml_map.at("pressure");
    pressure_r = primr_map.at("pressure");
  } else {
    pressure_l = EFlt3DArray(shape[0], shape[1], shape[2]);
    pressure_r = EFlt3DArray(shape[0], shape[1], shape[2]);
  }

}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) const throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Construct a map holding the field data for each of the (non-passive)
  // primitive quantities.
  EnzoEFltArrayMap primitive_map = nonpassive_primitive_map_(block);

  if (eos_->uses_dual_energy_formalism()){
    // synchronize eint and etot.
    // This is only strictly necessary after problem initialization and when
    // there is an inflow boundary condition
    eos_->apply_floor_to_energy_and_sync(primitive_map, 0);
  }

  // Compute thermal pressure (this presently requires that "pressure" is a
  // permanent field)
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray pressure = array_factory.from_name("pressure");
  EnzoEFltArrayMap conserved_passive_scalar_map =
      conserved_passive_scalar_map_(block);
  eos_->pressure_from_integrable(primitive_map, pressure,
				 conserved_passive_scalar_map, 0);

  // Now load other necessary quantities
  enzo_float gamma = eos_->get_gamma();
  EFlt3DArray density = primitive_map.at("density");
  EFlt3DArray velocity_x = primitive_map.at("velocity_x");
  EFlt3DArray velocity_y = primitive_map.at("velocity_y");
  EFlt3DArray velocity_z = primitive_map.at("velocity_z");

  const bool mhd = (mhd_choice_ != bfield_choice::no_bfield);
  EFlt3DArray bfieldc_x, bfieldc_y, bfieldc_z;
  if (mhd) {
    bfieldc_x = primitive_map.at("bfield_x");
    bfieldc_y = primitive_map.at("bfield_y");
    bfieldc_z = primitive_map.at("bfield_z");
  }

  // widths of cells
  EnzoBlock * enzo_block = enzo::block(block);
  double dx = enzo_block->CellWidth[0];
  double dy = enzo_block->CellWidth[1];
  double dz = enzo_block->CellWidth[2];

  // initialize
  double dtBaryons = ENZO_HUGE_VAL;

  // timestep is the minimum of 0.5 * dr_i/(abs(v_i)+cfast) for all dimensions.
  // dr_i and v_i are the the width of the cell and velocity along dimension i.
  // cfast = fast magnetosonic speed (Convention is to use max value:
  // cfast = (va^2+cs^2)

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {
	enzo_float bmag_sq = 0.0;
        if (mhd){
          bmag_sq = (bfieldc_x(iz,iy,ix) * bfieldc_x(iz,iy,ix) +
		     bfieldc_y(iz,iy,ix) * bfieldc_y(iz,iy,ix) +
		     bfieldc_z(iz,iy,ix) * bfieldc_z(iz,iy,ix));
        }
	// Using "Rationalized" Gaussian units (where the magnetic permeability
	// is mu=1 and pressure = B^2/2)
	// To convert B to normal Gaussian units, multiply by sqrt(4*pi)
	enzo_float inv_dens= 1./density(iz,iy,ix);
	double cfast = (double) std::sqrt(gamma * pressure(iz,iy,ix) *
					  inv_dens + bmag_sq * inv_dens);

	dtBaryons = std::min(dtBaryons,
			     dx/(std::fabs((double) velocity_x(iz,iy,ix)) +
				 cfast));
	dtBaryons = std::min(dtBaryons,
			     dy/(std::fabs((double) velocity_y(iz,iy,ix)) +
				 cfast));
	dtBaryons = std::min(dtBaryons,
			     dz/(std::fabs((double) velocity_z(iz,iy,ix))
				 + cfast));
      }
    }
  }
  // Multiply resulting dt by CourantSafetyNumber (for extra safety!).
  // This should be less than 0.5 for standard algorithm
  dtBaryons *= courant_;

  return dtBaryons;
}
