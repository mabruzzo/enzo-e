// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Method class

#ifndef PROBLEM_METHOD_HPP
#define PROBLEM_METHOD_HPP

class Refresh;
class Schedule;
class PassiveSchedule;

class Method : public PUP::able 
{
  /// @class    Method
  /// @ingroup  Method
  /// @brief    [\ref Method] Interface to an application method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method (double courant = 1.0) throw();

  /// Destructor
  virtual ~Method() throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Method);
  
  Method (CkMigrateMessage *m)
    : PUP::able(m),
    passive_schedule_(NULL),
    courant_(1.0),
    ir_post_(-1),
    neighbor_type_(neighbor_leaf)
   
  { }
      
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw() = 0; 

  /// Return the name of this Method
  virtual std::string name () throw () = 0;

  /// Compute maximum timestep for this method
  virtual double timestep (Block * block) const throw() 
  { return std::numeric_limits<double>::max(); }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw()
  {
    /* This function intentionally empty */
  }

  /// Add a new refresh object
  int add_new_refresh_ (int neighbor_type = neighbor_leaf);

  /// Return the index for the main post-refresh object
  int refresh_id_post() const;

  /// advance schedule (if necessary) and return whether method is scheduled
  bool advance_schedule(int cycle, double time) throw();

  /// Set schedule
  void set_schedule (Schedule * schedule, double initial_time) throw();

  double courant() const throw ()
  { return courant_; }

  void set_courant(double courant) throw ()
  { courant_ = courant; }

protected: // functions

  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y,
	      int mx, int my, int mz,
	      bool active = true) const throw()
  {
    if (! active ) return;
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) X[i] = Y[i];
  }

public: // attributes (static)

  static double courant_global;

protected: // attributes

  /// PassiveSchedule object, if any (default is every cycle)
  PassiveSchedule * passive_schedule_;

  /// Courant condition for the Method
  double courant_;

  /// Index for main refresh after Method is called
  int ir_post_;

  /// Default refresh type
  int neighbor_type_;
};

#endif /* PROBLEM_METHOD_HPP */
