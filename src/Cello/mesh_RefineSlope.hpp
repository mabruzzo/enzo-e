// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the RefineSlope class
///

#ifndef MESH_REFINE_SLOPE_HPP
#define MESH_REFINE_SLOPE_HPP

class RefineSlope : public Refine {
  /// @class    RefineSlope
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh]

public:  // interface
  /// Constructor
  RefineSlope(double min_refine, double max_coarsen,
              std::vector<std::string> field_name_list, int max_level,
              bool include_ghosts, std::string output) throw();

  /// default constructor
  // RefineSlope () throw() : Refine() {};

  PUPable_decl(RefineSlope);

  RefineSlope(CkMigrateMessage* m) : Refine(m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup(PUP::er& p) {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | field_id_list_;
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply(Block* block) throw();

  virtual std::string name() const { return "slope"; };

private:  // functions
  template <class T>
  void evaluate_block_(T* array, T* output, int ndx, int ndy, int ndz, int gx,
                       int gy, int gz, bool* any_refine, bool* all_coarsen,
                       int rank, double* h3);

private:  // attributes
  /// List of field id's
  std::vector<int> field_id_list_;
};

#endif /* MESH_REFINE_SLOPE_HPP */
