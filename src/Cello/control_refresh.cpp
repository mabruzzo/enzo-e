// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_REFRESH

#ifdef DEBUG_REFRESH
#  define TRACE_REFRESH(msg,REFRESH)				\
  printf ("%d %s:%d %s TRACE_REFRESH %s type %d\n",CkMyPe(),	\
	  __FILE__,__LINE__,name().c_str(),msg,REFRESH->sync_type());	\
  fflush(stdout);
#else
#  define TRACE_REFRESH(msg,REFRESH) /* NOTHING */
#endif

//----------------------------------------------------------------------

void Block::refresh_enter (int callback, Refresh * refresh) 
{
#ifdef DEBUG_REFRESH
  CkPrintf ("%d %s:%d DEBUG REFRESH %s Block::refresh_enter(%p) callback %d refresh callback %d\n",
	    CkMyPe(), __FILE__,__LINE__,name().c_str(),refresh,callback,refresh->callback());
  fflush(stdout);
#endif
  
  set_refresh(refresh);

  // Update refresh object for the Block

  refresh_.back()->set_callback(callback);

  refresh_begin_();
}

//----------------------------------------------------------------------

void Block::refresh_begin_() 
{

  Refresh * refresh = this->refresh();
  TRACE_REFRESH("refresh_begin_()",refresh);

  check_leaf_();

  check_delete_();

  cello::simulation()->set_phase(phase_refresh);

  control_sync (CkIndex_Block::p_refresh_continue(),
		refresh->sync_type(),
		refresh->sync_load(),
		refresh->min_face_rank(),
		refresh->neighbor_type(),
		refresh->root_level());
}

//----------------------------------------------------------------------

void Block::refresh_continue()
{

  // Refresh if Refresh object exists and have data

  Refresh * refresh = this->refresh();
  TRACE_REFRESH("refresh_continue_()",refresh);

  if ( refresh && refresh->active() ) {

    // count self
    int count = 1;

    // send Field face data
    if (refresh->any_fields()) {
      count += refresh_load_field_faces_ (refresh);
    }

    // send Particle face data
    if (refresh->any_particles()){
      count += refresh_load_particle_faces_ (refresh);
    }

    // wait for all messages to arrive (including this one)
    // before continuing to p_refresh_exit()

    control_sync_count(CkIndex_Block::p_refresh_exit(),
			refresh->sync_store(),count);
    
  } else {

    refresh_exit_();
    
  }
  
}

//----------------------------------------------------------------------


void Block::p_refresh_store (MsgRefresh * msg)
{

  
  performance_start_(perf_refresh_store);

  msg->update(data());

  delete msg;

  Refresh * refresh = this->refresh();
  TRACE_REFRESH("p_refresh_store()",refresh);

  control_sync_count(CkIndex_Block::p_refresh_exit(),
		     refresh->sync_store(),0);
  
  performance_stop_(perf_refresh_store);
  performance_start_(perf_refresh_store_sync);
}


//----------------------------------------------------------------------

int Block::refresh_load_field_faces_ (Refresh *refresh)
{

  int count = 0;

  const int min_face_rank = refresh->min_face_rank();
  const int neighbor_type = refresh->neighbor_type();

  if (neighbor_type == neighbor_leaf ||
      neighbor_type == neighbor_tree) {

    // Loop over neighbor leaf Blocks (not necessarily same level)

    const int min_level = cello::config()->mesh_min_level;
    
    ItNeighbor it_neighbor =
      this->it_neighbor(min_face_rank,index_,
			neighbor_type,min_level,refresh->root_level());

    int if3[3];
    while (it_neighbor.next(if3)) {

      Index index_neighbor = it_neighbor.index();

      int ic3[3];
      it_neighbor.child(ic3);

      const int level = this->level();
      const int level_face = it_neighbor.face_level();

      const int refresh_type = 
	(level_face == level - 1) ? refresh_coarse :
	(level_face == level)     ? refresh_same :
	(level_face == level + 1) ? refresh_fine : refresh_unknown;

      refresh_load_field_face_ (refresh_type,index_neighbor,if3,ic3);
      ++count;
    }

  } else if (neighbor_type == neighbor_level) {

    // Loop over neighbor Blocks in same level (not necessarily leaves)

    ItFace it_face = this->it_face(min_face_rank,index_);

    int if3[3];
    while (it_face.next(if3)) {

      // count all faces if not a leaf, else don't count if face level
      // is less than this block's level
      
      if ( ! is_leaf() || face_level(if3) >= level()) {
	
	Index index_face = it_face.index();
	int ic3[3] = {0,0,0};
	refresh_load_field_face_ (refresh_same,index_face,if3,ic3);
	++count;
      }

    }
  }

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_load_field_face_
( int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  //  TRACE_REFRESH("refresh_load_field_face()");

  // REFRESH FIELDS

  // ... coarse neighbor requires child index of self in parent

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object
  bool lg3[3] = {false,false,false};

  Refresh * refresh = this->refresh();

  FieldFace * field_face = create_face
    (if3, ic3, lg3, refresh_type, refresh,false);
#ifdef DEBUG_FIELD_FACE  
  CkPrintf ("%d %s:%d DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,__LINE__,field_face);
#endif

  DataMsg * data_msg = new DataMsg;

  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (data()->field_data(),false);

  MsgRefresh * msg = new MsgRefresh;

  msg->set_data_msg (data_msg);

  thisProxy[index_neighbor].p_refresh_store (msg);

}


//----------------------------------------------------------------------

int Block::refresh_load_particle_faces_ (Refresh * refresh)
{

  //  TRACE_REFRESH("refresh_load_particle_faces()");
  
  const int rank = cello::rank();

  const int npa = (rank == 1) ? 4 : ((rank == 2) ? 4*4 : 4*4*4);

  ParticleData * particle_array[npa];
  ParticleData * particle_list [npa];
  Index * index_list = new Index[npa];
  
  for (int i=0; i<npa; i++) {
    particle_list[i]  = NULL;
    particle_array[i] = NULL;
  }

  // Sort particles that have left the Block into 4x4x4 array
  // corresponding to neighbors

  int nl = particle_load_faces_
    (npa,particle_list,particle_array, index_list, refresh);

  // Send particle data to neighbors

  particle_send_(nl,index_list,particle_list);

  delete [] index_list;

  return nl;
}
//----------------------------------------------------------------------

int Block::particle_load_faces_ (int npa, 
				 ParticleData * particle_list[],
				 ParticleData * particle_array[],
				 Index index_list[],
				 Refresh *refresh)
{
  // Array elements correspond to child-sized blocks to
  // the left, inside, and right of the main Block.  Particles
  // are assumed to be (well) within this area.
  //
  //     +---+---+---+---+
  //     | 03| 13| 23| 33|
  //     +---+===+===+---+
  //     | 02||  :  || 32|
  //     +---+ - + - +---+
  //     | 01||  :  || 31|
  //     +---+=======+---+
  //     | 00| 10| 20| 30|
  //     +---+---+---+---+
  //
  // Actual neighbors may overlap multiple child-sized blocks.  In
  // that case, we have one ParticleData object per neighbor, but
  // with pointer duplicated.   So if neighbor configuration is:
  //
  //     +---+   5   +---+
  //     | 4 |       | 6 |
  // +---+---+===+===+---+
  // |       ||     ||    
  // |   2   +       +   3
  // |       ||     ||    
  // +-------+=======+-------+
  //         |            
  //     0   |            
  //                 1   
  //
  // Then the particle data array will be:
  //
  //     +---+---+---+---+
  //     | 4 | 5 | 5 | 6 |
  //     +---+===+===+---+
  //     | 2 ||  :  || 3 |
  //     +---+ - + - +---+
  //     | 2 ||  :  || 3 |
  //     +---+=======+---+
  //     | 0 | 1 | 1 | 1 |
  //     +---+---+---+---+

  //  TRACE_REFRESH("particle_load_faces()");

  // ... arrays for updating positions of particles that cross
  // periodic boundaries

  int nl = particle_create_array_neighbors_
    (refresh, particle_array,particle_list,index_list);

  // Scatter particles among particle_data array

  Particle particle (cello::particle_descr(),
		     data()->particle_data());

  std::vector<int> type_list;
  if (refresh->all_particles()) {
    const int nt = particle.num_types();
    type_list.resize(nt);
    for (int i=0; i<nt; i++) type_list[i] = i;
  } else {
    type_list = refresh->particle_list();
  }

  particle_scatter_neighbors_(npa,particle_array,type_list, particle);

  // Update positions particles crossing periodic boundaries

  particle_apply_periodic_update_  (nl,particle_list,refresh);

  return nl;
}

//----------------------------------------------------------------------

int Block::particle_create_array_neighbors_
(Refresh * refresh, 
 ParticleData * particle_array[],
 ParticleData * particle_list[],
 Index index_list[])
{ 
  //  TRACE_REFRESH("particle_create_array_neighbors()");

  const int rank = cello::rank();
  const int level = this->level();

  const int min_face_rank = refresh->min_face_rank();
  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_,neighbor_leaf,0,0);

  int il = 0;

  int if3[3];
  for (il=0; it_neighbor.next(if3); il++) {

    const int level_face = it_neighbor.face_level();

    int ic3[3] = {0,0,0};

    const int refresh_type = 
      (level_face == level - 1) ? refresh_coarse :
      (level_face == level)     ? refresh_same :
      (level_face == level + 1) ? refresh_fine : refresh_unknown;

    if (refresh_type==refresh_coarse) {
      // coarse neighbor: need index of self in parent
      index_.child(index_.level(),ic3,ic3+1,ic3+2);
    } else if (refresh_type==refresh_fine) {
      // fine neighbor: need index of child in self
      it_neighbor.child(ic3);
    }
    // (else same-level neighbor: don't need child)

    int index_lower[3] = {0,0,0};
    int index_upper[3] = {1,1,1};
    refresh->index_limits (rank,refresh_type,if3,ic3,index_lower,index_upper);

    ParticleData * pd = new ParticleData;

    ParticleDescr * p_descr = cello::particle_descr();

    pd->allocate(p_descr);

    particle_list[il] = pd;

    index_list[il] = it_neighbor.index();

    for (int iz=index_lower[2]; iz<index_upper[2]; iz++) {
      for (int iy=index_lower[1]; iy<index_upper[1]; iy++) {
	for (int ix=index_lower[0]; ix<index_upper[0]; ix++) {
	  int i=ix + 4*(iy + 4*iz);
	  particle_array[i] = pd;
	}
      }
    }
  }
  
  return il;
}

//----------------------------------------------------------------------

void Block::particle_determine_periodic_update_
(int * index_lower, int * index_upper,
 double *dpx, double *dpy, double *dpz)
{
  //     ... domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;

  cello::hierarchy()->lower(&dxm,&dym,&dzm);
  cello::hierarchy()->upper(&dxp,&dyp,&dzp);

  //     ... periodicity
  bool p3[3];
  periodicity(p3);

  //     ... boundary
  bool b32[3][2];
  is_on_boundary (b32);

  const int rank = cello::rank();

  // Update (dpx,dpy,dpz) position correction if periodic domain
  // boundary is crossed

  if (rank >= 1) {
    if (index_lower[0]==0 && b32[0][0] && p3[0]) (*dpx) = +(dxp - dxm);
    if (index_upper[0]==4 && b32[0][1] && p3[0]) (*dpx) = -(dxp - dxm);
  }
  if (rank >= 2) {
    if (index_lower[1]==0 && b32[1][0] && p3[1]) (*dpy) = +(dyp - dym);
    if (index_upper[1]==4 && b32[1][1] && p3[1]) (*dpy) = -(dyp - dym);
  }
  if (rank >= 3) {
    if (index_lower[2]==0 && b32[2][0] && p3[2]) (*dpz) = +(dzp - dzm);
    if (index_upper[2]==4 && b32[2][1] && p3[2]) (*dpz) = -(dzp - dzm);
  }
}

//----------------------------------------------------------------------

void Block::particle_apply_periodic_update_
(int nl, ParticleData * particle_list[], Refresh * refresh)
{

  const int rank = cello::rank();
  const int level = this->level();
  const int min_face_rank = refresh->min_face_rank();

  double dpx[nl],dpy[nl],dpz[nl];

  for (int i=0; i<nl; i++) {
    dpx[i]=0.0;
    dpy[i]=0.0;
    dpz[i]=0.0; 
  }

  // Compute position updates for particles crossing periodic boundaries

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_,neighbor_leaf,0,0);

  int il=0;

  int if3[3];
  while (it_neighbor.next(if3)) {

    const int level_face = it_neighbor.face_level();

    int ic3[3];
    it_neighbor.child(ic3);

    const int refresh_type = 
      (level_face == level - 1) ? refresh_coarse :
      (level_face == level)     ? refresh_same :
      (level_face == level + 1) ? refresh_fine : refresh_unknown;

    int index_lower[3] = {0,0,0};
    int index_upper[3] = {1,1,1};
    refresh->index_limits (rank,refresh_type,if3,ic3,index_lower,index_upper);

    // ASSERT: il < nl
    particle_determine_periodic_update_
      (index_lower,index_upper,&dpx[il],&dpy[il],&dpz[il]);

    il++;

  }

  ParticleDescr * p_descr = cello::particle_descr();

  // Apply the updates to the list of particles

  for (int il=0; il<nl; il++) {

    ParticleData * p_data = particle_list[il];
    Particle particle_neighbor (p_descr,p_data);

    if ( ((rank >= 1) && dpx[il] != 0.0) ||
	 ((rank >= 2) && dpy[il] != 0.0) ||
	 ((rank >= 3) && dpz[il] != 0.0) ) {
	
      // ... for each particle type
      const int nt = particle_neighbor.num_types();
      for (int it=0; it<nt; it++) {

	// ... for each batch of particles
	const int nb = particle_neighbor.num_batches(it);
	for (int ib=0; ib<nb; ib++) {

	  particle_neighbor.position_update (it,ib,dpx[il],dpy[il],dpz[il]);

	}
      }
    }
  }
}
//----------------------------------------------------------------------

void Block::particle_scatter_neighbors_
(int npa,
 ParticleData * particle_array[],
 std::vector<int> & type_list,
 Particle particle)
{
  const int rank = cello::rank();

  //     ... get Block bounds 
  double xm,ym,zm;
  double xp,yp,zp;
  lower(&xm,&ym,&zm);
  upper(&xp,&yp,&zp);

  // find block center (x0,y0,z0) and width (xl,yl,zl)
  const double x0 = 0.5*(xm+xp);
  const double y0 = 0.5*(ym+yp);
  const double z0 = 0.5*(zm+zp);
  const double xl = xp-xm;
  const double yl = yp-ym;
  const double zl = zp-zm;

  int count = 0;
  // ...for each particle type to be moved

  for (auto it_type=type_list.begin(); it_type!=type_list.end(); it_type++) {

    int it = *it_type;

    const int ia_x  = particle.attribute_position(it,0);

    // (...positions may use absolute coordinates (float) or
    // block-local coordinates (int))
    const bool is_float = 
      (cello::type_is_float(particle.attribute_type(it,ia_x)));

    // (...stride may be != 1 if particle attributes are interleaved)
    const int d  = particle.stride(it,ia_x);

    // ...for each batch of particles

    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);

      // ...extract particle position arrays

      double xa[np],ya[np],za[np];
      particle.position(it,ib,xa,ya,za);

      // ...initialize mask used for scatter and delete
      // ...and corresponding particle indices

      bool mask[np];
      int index[np];
      for (int ip=0; ip<np; ip++) {

	double x = is_float ? 2.0*(xa[ip*d]-x0)/xl : xa[ip*d];
	double y = is_float ? 2.0*(ya[ip*d]-y0)/yl : ya[ip*d];
	double z = is_float ? 2.0*(za[ip*d]-z0)/zl : za[ip*d];

	int ix = (rank >= 1) ? (x + 2) : 0;
	int iy = (rank >= 2) ? (y + 2) : 0;
	int iz = (rank >= 3) ? (z + 2) : 0;

	if (! (0 <= ix && ix < 4) ||
	    ! (0 <= iy && iy < 4) ||
	    ! (0 <= iz && iz < 4)) {
	  
	  CkPrintf ("%d ix iy iz %d %d %d\n",CkMyPe(),ix,iy,iz);
	  CkPrintf ("%d x y z %f %f %f\n",CkMyPe(),x,y,z);
	  CkPrintf ("%d xa ya za %f %f %f\n",CkMyPe(),xa[ip*d],ya[ip*d],za[ip*d]);
	  CkPrintf ("%d xm ym zm %f %f %f\n",CkMyPe(),xm,ym,zm);
	  CkPrintf ("%d xp yp zp %f %f %f\n",CkMyPe(),xp,yp,zp);
	  ERROR3 ("Block::particle_scatter_neighbors_",
		  "particle indices (ix,iy,iz) = (%d,%d,%d) out of bounds",
		  ix,iy,iz);
	}

	const int i = ix + 4*(iy + 4*iz);
	index[ip] = i;
	bool in_block = true;
	in_block = in_block && (!(rank >= 1) || (1 <= ix && ix <= 2));
	in_block = in_block && (!(rank >= 2) || (1 <= iy && iy <= 2));
	in_block = in_block && (!(rank >= 3) || (1 <= iz && iz <= 2));
	mask[ip] = ! in_block;
      }

      // ...scatter particles to particle array
      particle.scatter (it,ib, np, mask, index, npa, particle_array);
      // ... delete scattered particles
      count += particle.delete_particles (it,ib,mask);
    }
  }

  cello::simulation()->data_delete_particles(count);

}

//----------------------------------------------------------------------

void Block::particle_send_
(int nl,Index index_list[], ParticleData * particle_list[])
{

  ParticleDescr * p_descr = cello::particle_descr();

  for (int il=0; il<nl; il++) {

    Index index           = index_list[il];
    ParticleData * p_data = particle_list[il];
    Particle particle_send (p_descr,p_data);
    
    if (p_data && p_data->num_particles(p_descr)>0) {

      DataMsg * data_msg = new DataMsg;
      data_msg ->set_particle_data(p_data,true);

      MsgRefresh * msg = new MsgRefresh;
      msg->set_data_msg (data_msg);

      thisProxy[index].p_refresh_store (msg);

    } else if (p_data) {
      
      MsgRefresh * msg = new MsgRefresh;

      thisProxy[index].p_refresh_store (msg);

      // assert ParticleData object exits but has no particles
      delete p_data;

    }

  }
}
