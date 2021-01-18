// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataMsg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-06
/// @brief    [\ref Data] Declaration of the DataMsg Charm++ Message

#ifndef DATA_DATA_MSG_HPP
#define DATA_DATA_MSG_HPP

class FaceFluxes;
class ParticleData;
class FieldData;
class FieldFace;

class DataMsg {

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  DataMsg() 
    : field_face_(nullptr),
      field_face_delete_   (false),
      field_data_(nullptr),
      field_data_delete_   (false),
      particle_data_(nullptr),
      particle_data_delete_(false),
      face_fluxes_list_(),
      face_fluxes_delete_(),
      padded_face_(),
      padded_face_field_list_()
  {
    cello::hex_string(tag_,8);
    for (int i=0; i<3; i++) {
      ma3_pf_[i] =0;
      n3_pf_[i]  =0;
      if3_pf_[i] =0;
      iam3_pf_[i]=0;
    }
    ++counter[cello::index_static()]; 
  }

  ~DataMsg()
  {
    --counter[cello::index_static()];
    
    if (field_face_delete_) {
      delete field_face_;
      field_face_ = nullptr;
    }
    if (field_data_delete_) {
      delete field_data_;
      field_data_ = nullptr;
    }
    if (particle_data_delete_) {
      delete particle_data_;
      particle_data_ = nullptr;
    }
    for (size_t i=0; i<face_fluxes_list_.size(); i++) {
      if (face_fluxes_delete_[i]) {
        delete face_fluxes_list_[i];
        face_fluxes_list_[i] = nullptr;
      }
    }
    face_fluxes_list_.clear();
    face_fluxes_delete_.clear();
    padded_face_.clear();
    padded_face_field_list_.clear();
  }

  /// Copy constructor
  DataMsg(const DataMsg & data_msg) throw()
  {
    ++counter[cello::index_static()]; 
  };

  /// Assignment operator
  DataMsg & operator= (const DataMsg & data_msg) throw()
  { return *this; }

  void pup(PUP::er &p) {
    TRACEPUP;
    WARNING("DataMsg::pup()",
	    "DataMsg::pup() should never be called!");
  }

  /// --------------------
  /// FIELD DATA
  /// --------------------

  /// Return the FieldFace
  FieldFace * field_face () 
  {
    return field_face_;
  }

  /// Set the FieldFace object
  void set_field_face  (FieldFace * field_face, bool is_new) 
  {
    field_face_ = field_face; 
    field_face_delete_ = is_new;
  }

  /// Return the serialized FieldFace array
  char * field_array () 
  { return field_array_; }


  /// Return the FieldData
  FieldData * field_data () 
  { return field_data_; }

  /// Set the FieldData object
  void set_field_data    (FieldData * field_data, bool is_new) 
  {
    field_data_ = field_data;
    field_data_delete_ = is_new;
  }

  /// --------------------
  /// PARTICLE DATA
  /// --------------------
  
  /// Return the ParticleData
  ParticleData * particle_data () 
  { return particle_data_; }

  /// Set the ParticleData object
  void set_particle_data  (ParticleData * particle_data, bool is_new) 
  {
    particle_data_ = particle_data; 
    particle_data_delete_ = is_new;
  }

  /// Delete the ParticleData object
  void delete_particle_data  () 
  { 
    delete particle_data_;
    particle_data_ = nullptr; 
  }

  /// --------------------
  /// FLUX DATA
  /// --------------------
  
  /// Return the ith FaceFluxes
  FaceFluxes * face_fluxes (unsigned i) 
  { return face_fluxes_list_[i]; }

  /// Return the number of FaceFluxes
  unsigned num_face_fluxes() const
  { return face_fluxes_list_.size(); }

  void set_num_face_fluxes(unsigned i)
  {
    if (i > face_fluxes_list_.size()) {
      face_fluxes_list_.resize(i);
      face_fluxes_delete_.resize(i);
    }
  }
  
  /// Set the FaceFluxes object
  void set_face_fluxes  (unsigned i, FaceFluxes * face_fluxes, bool is_new) 
  {
    set_num_face_fluxes(i+1);
    face_fluxes_list_[i] = face_fluxes; 
    face_fluxes_delete_[i] = is_new;
  }

  /// ------------------
  /// PADDED FACE ARRAYS
  /// ------------------

  /// Initialize the padded face arrays to send to neighbors
  void set_padded_face
  (int ma3[3], int n3[3], int r, int v,
   int iam3[3], int ifm3[3], int if3[3],
   std::vector<int> padded_face_field_list_, Field field,
   std::string block_name_temp);

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------
  
  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

  /// Update the Data with the data stored in this DataMsg
  void update (Data * data, bool is_local);

  /// Debugging
  void print (const char * message) const;

public: // static methods

  
protected: // attributes

  /// Field Face Data
  FieldFace * field_face_;
  /// Whethere FieldFace data should be deleted in destructor
  bool field_face_delete_;

  /// Field data
  union {

    /// Field data if local
    FieldData * field_data_;

    /// packed source field data if remote
    char * field_array_;

  };
  /// Whethere FieldData data should be deleted in destructor
  bool field_data_delete_;
  
  /// Particle data
  ParticleData * particle_data_;
  /// Whethere Particle data should be deleted in destructor
  bool particle_data_delete_;

  /// Flux faces (array for each field)
  std::vector<FaceFluxes *> face_fluxes_list_;

  /// Whether Flux data should be deleted in destructor
  std::vector<bool> face_fluxes_delete_;

  /// Padded face array (for interpolation that requires extra layer
  /// of cells around FieldFace that intersects multiple Blocks)
  /// This stores padded faces for all fields in
  std::vector<cello_float> padded_face_;
  /// List of field indices for padded face
  std::vector<int> padded_face_field_list_;
  /// dimensions of the padded face array
  int ma3_pf_[3];
  /// size of the coarse-block array section
  int n3_pf_[3];
  /// lower loop limits for the padded array
  int iam3_pf_[3];
  /// face associated with the padded face (from the receiver's perspective)
  int if3_pf_[3];
  /// hex tag identifying object to match sender and receiver
  char tag_[9];

};

#endif /* DATA_DATA_MSG_HPP */

