// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-09-04
/// @brief

#include "problem.hpp"

double Method::courant_global = 1.0;

//----------------------------------------------------------------------

Method::Method (double courant) throw()
  : schedule_(NULL),
    courant_(courant),
    neighbor_type_(neighbor_leaf)
{
  ir_post_ = add_new_refresh_();
  cello::refresh(ir_post_)->set_callback(CkIndex_Block::p_compute_continue());
}

//----------------------------------------------------------------------

Method::~Method() throw()
{
  delete schedule_;
}

//----------------------------------------------------------------------

void Method::pup (PUP::er &p)
{ TRACEPUP;
  PUP::able::pup(p);

  p | schedule_; // pupable
  p | courant_;
  p | ir_post_;
  p | neighbor_type_;
  p | required_fields_; // std::vector<str> required fields
  p | field_centering_; // std::map<std::string, std::array<int,3>>

}

//----------------------------------------------------------------------

int Method::add_new_refresh_ (int neighbor_type)
{
  // set Method::ir_post_

  const int ghost_depth = 4; // std::max(g3[0],std::max(g3[1],g3[2]));
  const int min_face_rank = 0; // cello::config()->adapt_min_face_rank;

  // Set default refresh object
  Refresh refresh_default
    (ghost_depth,min_face_rank, neighbor_type, sync_neighbor, 0);

  return cello::simulation()->new_register_refresh(refresh_default);
}

//----------------------------------------------------------------------

int Method::refresh_id_post() const
{
  return ir_post_;
}

//----------------------------------------------------------------------

void Method::set_schedule (Schedule * schedule) throw()
{
  if (schedule_) delete schedule_;
  schedule_ = schedule;
}

//----------------------------------------------------------------------

void Method::define_fields () throw()
{
  /* Ensure required fields are defined for this method */

  FieldDescr * field_descr = cello::field_descr();
  Config   * config  = (Config *) cello::config();

  bool added_fields = false;

  for (const std::string &field : required_fields_){

    // get the field_centering
    int cx, cy, cz;
    
    if ( field_centering_.find(field) != field_centering_.end()){
      cx = field_centering_[field][0];
      cy = field_centering_[field][1];
      cz = field_centering_[field][2];
    } else { // the field is assumed to be cell-centered, by default
      cx = 0;
      cy = 0;
      cz = 0;
    }

    if( ! field_descr->is_field( field )){ // insert new field
      int id_field = field_descr->insert_permanent( field );
      field_descr->set_precision(id_field, config->field_precision);
      field_descr->set_centering(id_field, cx, cy, cz);

      added_fields = true;

    } else { // validate the centering of the existing field
      int id_field = field_descr->field_id(field);
      int cur_cx, cur_cy, cur_cz;
      field_descr->centering(id_field, &cur_cx, &cur_cy, &cur_cz);

      ASSERT7("Method::define_fields",
              ("The \"%s\" field is required to have a centering (cx,cy,cz) "
               "of (%d,%d,%d). It's pre-defined with a centering (%d,%d,%d)"),
              field.c_str(),   cx,cy,cz,   cur_cx,cur_cy,cur_cz,
              (cx==cur_cx) && (cy==cur_cy) && (cz==cur_cz));
    }
  }

  // Need to reconstruct history if new fields added
  if (added_fields) field_descr->reset_history(config->field_history);
}

//----------------------------------------------------------------------

void Method::define_group_fields (std::vector<std::string> group_fields,
                                  std::string groupname) throw()
{
  /* Ensure fields are grouped correctly */

  FieldDescr * field_descr = cello::field_descr();
  Config   * config  = (Config *) cello::config();

  bool added_fields = false;

  for (int ifield = 0; ifield < group_fields.size(); ifield++){

    // Maybe just throw error here to keep this fully separate from above
    if( ! field_descr->is_field( required_fields_[ifield] )){
      int field_id = field_descr->insert_permanent( required_fields_[ifield] );
      field_descr->set_precision(field_id, config->field_precision);
      added_fields = true;
    }

    if (!(field_descr->groups()->is_in( group_fields[ifield], groupname)) ){
      field_descr->groups()->add( group_fields[ifield], groupname);
    }

  }

  // Need to reconstruct history if new fields added
  if (added_fields) field_descr->reset_history(config->field_history);
}

//======================================================================
