// See LICENSE_CELLO file for license and copyright information

/// @file     memory_Memory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Sep  3 16:29:56 PDT 2009
/// @brief    [\ref Memory] Interface for the Memory class.  Uses the
/// Singleton design pattern.

#ifndef MEMORY_MEMORY_HPP
#define MEMORY_MEMORY_HPP

class Memory {
  /// @class    Memory
  /// @ingroup  Memory
  /// @brief    [\ref Memory] Manage memory allocation and deallocation

public:  // interface
         /// Get single instance of the Memory object
  static Memory* instance() {
#ifdef CONFIG_USE_MEMORY
    return &instance_[cello::index_static()];
#else
    return 0;
#endif
  }

private:  // interface
  /// Create the (single) Memory object (singleton design pattern)
  Memory()
#ifdef CONFIG_USE_MEMORY
      : is_active_(false),
        warning_mb_(0.0),
        limit_gb_(0.0)
#endif
  {
    initialize_();
  }

  /// Copy the (single) Memory object (singleton design pattern)
  Memory(const Memory&);

  /// Assign the (single) Memory object (singleton design pattern)
  Memory& operator=(const Memory& memory);

  /// Delete the Memory object
  ~Memory() {}

public:  // interface
  void pup(PUP::er& p) {
    TRACEPUP;
#ifdef CONFIG_USE_MEMORY
    p | is_active_;
    p | fill_new_;
    p | fill_delete_;
    p | bytes_limit_;
#endif
    static bool warn[CONFIG_NODE_SIZE] = {false};
    const int in = cello::index_static();
    if (!warn[in]) {
      WARNING("Memory::pup()", "Skipping index_group_");
      warn[in] = true;
    }
    p | group_name_;
    if (p.isUnpacking()) {
      WARNING("Memory::pup()", "Calling Memory::reset()");
      reset();
    }
  }

  /// Allocate memory
  void* allocate(size_t size);

  /// De-allocate memory
  void deallocate(void* pointer);

  /// Define a new group
  void new_group(std::string group_name);

  int index_group(std::string group_name = "") const {
    int index_group = 0;
    for (size_t index = 0; index < group_name_.size(); index++) {
      if (group_name_.at(index) == group_name) index_group = index;
    }
    return index_group;
  }

  /// Begin allocating/deallocating memory associated with the given group
  void set_group(std::string group_name) {
    index_group_ = this->index_group(group_name);
  }

  std::string group() const { return group_name_[index_group_]; }

  /// Current number of bytes allocated
  int64_t bytes(std::string group = "");

  /// Estimate of used / available memory
  float efficiency(std::string group = "");

  /// Maximum number of bytes allocated within an interval
  int64_t bytes_high(std::string group = "");

  /// Maximum number of bytes allocated during run
  int64_t bytes_highest(std::string group = "");

  /// Specify the maximum number of bytes to use
  void set_bytes_limit(int64_t size, std::string group = "");

  /// Return the limit on number of bytes to use
  int64_t bytes_limit(std::string group = "");

  /// Query the maximum number of bytes left to use, 0 if unknown
  int64_t bytes_available(std::string group = "");

  /// Return the number of calls to allocate for the group
  int num_new(std::string group = "");

  /// Return the number of calls to deallocate for the group
  int num_delete(std::string group = "");

  /// Print memory summary
  void print();

  /// Reset memory counters for the current group
  void reset();

  /// Reset bytes_high to current
  void reset_high();

  /// Set whether memory tracking is active or not
  void set_active(bool is_active)
#ifdef CONFIG_USE_MEMORY
  {
    is_active_ = is_active;
  }
#else
  {
  }
#endif

  bool is_active() const {
#ifdef CONFIG_USE_MEMORY
    return is_active_;
#else
    return 0;
#endif
  }

  /// Set whether to fill memory with a value after allocating
  void set_fill_new(unsigned char value = 0) {
#ifdef CONFIG_USE_MEMORY
    fill_new_ = value;
#endif
  }

  /// Set whether to fill memory with a value before deallocating
  void set_fill_delete(unsigned char value = 0) {
#ifdef CONFIG_USE_MEMORY
    fill_delete_ = value;
#endif
  }

  void set_warning_mb(float value) {
#ifdef CONFIG_USE_MEMORY
    warning_mb_ = value;
#endif
  }

  void set_limit_gb(float value) {
#ifdef CONFIG_USE_MEMORY
    limit_gb_ = value;
#endif
  }

  //======================================================================

private:  // functions
  /// Initialize the memory component
  void initialize_();

  //======================================================================

private:  // attributes
#ifdef CONFIG_USE_MEMORY
          /// Single instance of the Memory object (singleton design pattern)
  static Memory instance_[CONFIG_NODE_SIZE];
#endif

#ifdef CONFIG_USE_MEMORY

  /// Whether keeping track of memory statistics is active or not
  bool is_active_;

  /// allocate clear value
  unsigned char fill_new_;

  /// deallocate clear value
  unsigned char fill_delete_;

  /// Limit on number of bytes to allocate
  std::vector<int64_t> bytes_limit_;

  /// Current bytes allocated for different groups
  std::vector<int64_t> bytes_curr_;

  /// Intervaled high-water bytes allocated for different groups
  std::vector<int64_t> bytes_high_;

  /// High-water bytes allocated for different groups
  std::vector<int64_t> bytes_highest_;

  /// Number of calls to new for different groups
  std::vector<int64_t> new_calls_;

  /// Number of calls to delete for different groups
  std::vector<int64_t> delete_calls_;

  /// Limit on single memory allocation to display warning
  float warning_mb_;

  /// Limit on total memory allocated before error (to prevent crashing machine)
  float limit_gb_;

#endif

  /// The current group index, or 0 if none
  int index_group_;

  /// Names of groups
  std::vector<std::string> group_name_;
};

#endif /* MEMORY_MEMORY_HPP */
