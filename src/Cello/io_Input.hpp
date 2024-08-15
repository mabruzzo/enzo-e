// See LICENSE_CELLO file for license and copyright information

/// @file     io_Input.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-23
/// @brief    [\ref Io] Declaration of the Input class

#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

class Factory;
class Hierarchy;
class ItIndex;
class Simulation;

class Input : public PUP::able

{
  /// @class    Input
  /// @ingroup  Io
  /// @brief [\ref Io] define interface for various types of IO for
  /// Simulations

public:  // functions
  /// empty constructor for charm++ pup()
  Input() throw() {}

  /// Create an uninitialized Input object
  Input(const Factory* factory) throw();

  /// Delete an Input object
  virtual ~Input() throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Input);

  /// Charm++ PUP::able migration constructor
  Input(CkMigrateMessage* m)
      : PUP::able(m),
        file_(0),
        sync_(0),
        index_charm_(0),
        cycle_(0),
        time_(0),
        file_name_(""),
        file_args_(),
        io_block_(0),
        it_field_index_(0),
        io_field_data_(0),
        it_particle_index_(0),
        io_particle_data_(0),
        stride_reader_(1) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p);

  /// Set file name
  void set_filename(std::string filename,
                    std::vector<std::string> fileargs) throw();

  /// Set field iterator
  void set_it_field_index(ItIndex* it_index) throw() {
    it_field_index_ = it_index;
  }

  /// Set particle iterator
  void set_it_particle_index(ItIndex* it_index) throw() {
    it_particle_index_ = it_index;
  }

  /// Return the IoBlock object
  IoBlock* io_block() const throw() { return io_block_; }

  /// Return the IoFieldData object
  IoFieldData* io_field_data() const throw() { return io_field_data_; }

  /// Return the IoParticleData object
  IoParticleData* io_particle_data() const throw() { return io_particle_data_; }

  /// Return the File object pointer
  File* file() throw() { return file_; };

  int stride_reader() const throw() { return stride_reader_; };

  void set_stride_reader(int stride) throw() {
    stride_reader_ = stride;
    sync_.set_stop(stride_reader_);
  };

  /// Accessor function for the CHARM Sync class
  Sync* sync() { return &sync_; };

  /// Set the index of this input in its simulation
  void set_index_charm(int index_charm) { index_charm_ = index_charm; }

public:  // virtual functions
  /// Initialize next input
  virtual void init() throw() {};

  /// Open (or create) a file for IO
  virtual void open() throw() = 0;

  /// Whether the file is open or not
  virtual bool is_open() throw() = 0;

  /// Close file for IO
  virtual void close() throw() = 0;

  /// Finalize input
  virtual void finalize() throw() {}

  /// Read metadata from the file
  void read_meta(Io* io) throw() { read_meta_(meta_type_file, io); }

  /// Read metadata from the current group in the file
  void read_meta_group(Io* io) throw() { read_meta_(meta_type_group, io); }

public:
  /// Read an entire simulation from disk
  virtual void read_simulation(Simulation* simulation) throw();

  /// Read local hierarchy data from disk
  virtual void read_hierarchy(Hierarchy* hierarchy) throw();

  /// Read local data data from disk
  virtual Block* read_block(Block* block, std::string block_name) throw();

  /// Read local field from disk
  virtual void read_field(Block* block, int index_field) throw() = 0;

  /// Read local particle from disk
  virtual void read_particle(Block* block, int index_particle) throw() = 0;

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote(int* n, char** buffer) throw() {};

  /// Accumulate and read data sent from a remote processes
  virtual void update_remote(int n, char* buffer) throw() {};

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote(int* n, char** buffer) throw() {};

protected:
  /// Return the filename for the file format and given arguments
  std::string expand_name_(const std::string* file_name,
                           const std::vector<std::string>* file_args) const
      throw();

private:
  void read_meta_(meta_type type, Io* io) throw();

protected:  // attributes
  /// File object for input
  File* file_;

  /// Sync for ending input
  Sync sync_;

  /// Index of this Input object in Simulation
  size_t index_charm_;

  /// Simulation cycle for next IO
  int cycle_;

  /// Simulation time for next IO
  double time_;

  /// Name of the file to read, including format arguments
  std::string file_name_;

  /// Format strings for file name, if any ("cycle", "time", etc.)
  std::vector<std::string> file_args_;

  /// I/O Block data accessor
  IoBlock* io_block_;

  /// Iterator over field indices
  ItIndex* it_field_index_;

  /// I/O FieldData data accessor
  IoFieldData* io_field_data_;

  /// Iterator over particle indices
  ItIndex* it_particle_index_;

  /// I/O ParticleData data accessor
  IoParticleData* io_particle_data_;

  /// Only processes with id's divisible by stride_reader_ reads
  /// (1: all processes read; 2: 0,2,4,... read; np: root process reads)
  int stride_reader_;
};

#endif /* IO_INPUT_HPP */
