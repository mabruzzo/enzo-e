// See LICENSE_CELLO file for license and copyright information

/// @file     data_Grouping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 16 15:27:42 PDT 2014
/// @brief    [\ref Data] Declaration of the Grouping class
///
/// This class serves to define groups of Fields, Particles, etc. into
/// named categories.  For example, one can define a "color" group
/// for fields, and add all color fields to the "color" grouping.  The
/// API supports adding groups, adding "items" (field names, particle
/// set names, etc.) to groups, testing whether an item is included
/// in a group, returning the size of the group, and returning an
/// iterator for items in a group.
///
/// NOTE: This class is not named "Group" since it conflicts with a
/// Charm++ class by the same name

#ifndef DATA_GROUPING_HPP
#define DATA_GROUPING_HPP

class Grouping {

  /// @class    Grouping
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  // /// Constructor
  // Grouping() throw();

  // /// Copy constructor
  // Grouping(const Grouping & Grouping) throw();

  // /// Assignment operator
  // Grouping & operator= (const Grouping & Grouping) throw();

  // /// Destructor
  // ~Grouping() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { p | groups_; }

  //----------------------------------------------------------------------

  /// Add an item to a group
  void add(std::string item, std::string group) 
    throw()
  {
    std::pair<std::string,std::string> value(item,group);
    groups_.insert(value);
  }

  /// Return whether the item is in the given group
  bool is_in(std::string item, std::string group) const
    throw()
  {
    std::pair<std::string,std::string> value(item,group);
    return groups_.find(value) != groups_.end();
  }

  /// Return the number of items in the group
  int size(std::string item) const
  {
    int count = 0;
    for (auto it=groups_.begin(); it != groups_.end(); it++) {
      if (it->second == item) ++count;
    }
    return count;
  }

  /// Return the nth item in the group
  std::string item (std::string group, int index)
  {
    int count = 0;

    for (auto it=groups_.begin(); it != groups_.end(); it++) {
      if (it->second == group) {
	if (count == index) return it->first;
	++count;
      }
    }
    return "";
  }
protected: // functions

  // NOTE: change pup() function whenever attributes change

  std::set<std::pair<std::string,std::string> > groups_;

};

#endif /* DATA_GROUPING_HPP */

