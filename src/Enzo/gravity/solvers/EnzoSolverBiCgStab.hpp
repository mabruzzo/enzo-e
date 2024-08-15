// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverBiCgStab.hpp
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoSolverBiCgStab
///
/// Bicongugate gradient stabilized solver (BiCgStab) for solving
/// linear systems on field data.

#ifndef ENZO_ENZO_SOLVER_BICGSTAB_HPP
#define ENZO_ENZO_SOLVER_BICGSTAB_HPP

class EnzoSolverBiCgStab : public Solver {
  enum bcg {
    bcg_undefined,
    bcg_start_2,
    bcg_loop_0a,
    bcg_loop_6,
    bcg_loop_12,
    bcg_loop_14
  };

  /// @class    EnzoSolverBiCgStab
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] This class implements the BiCgStab Krylov
  /// linear solver.  Alone, this is more applicable to smaller
  /// problems since the solver doesn't scale as well as some other
  /// solvers (FFT, MG, etc.) for larger problems.  Alternately, a
  /// more scalable solver may be combined as a preconditioner for a
  /// robust and scalable overall solver.

public:  // interface
  /// normal constructor
  EnzoSolverBiCgStab(std::string name, std::string field_x, std::string field_b,
                     int monitor_iter, int restart_cycle, int solve_type,
                     int index_prolong, int index_restrict, int min_level,
                     int max_level, int iter_max, double res_tol,
                     int index_precon, int coarse_level);

  /// default constructor
  EnzoSolverBiCgStab()
      : Solver(),
        A_(nullptr),
        is_alpha_(-1),
        is_beta_n_(-1),
        is_beta_d_(-1),
        is_rho0_(-1),
        is_err_(-1),
        is_err0_(-1),
        is_err_min_(-1),
        is_err_max_(-1),
        is_omega_(-1),
        is_omega_n_(-1),
        is_omega_d_(-1),
        is_rr_(-1),
        is_r0s_(-1),
        is_c_(-1),
        is_bs_(-1),
        is_xs_(-1),
        is_bnorm_(-1),
        is_vr0_(-1),
        is_ys_(-1),
        is_vs_(-1),
        is_us_(-1),
        is_qs_(-1),
        is_dot_sync_(-1),
        is_iter_(-1),
        res_tol_(0),
        index_precon_(-1),
        iter_max_(-1),
        ir_(-1),
        ir0_(-1),
        ip_(-1),
        iy_(-1),
        iv_(-1),
        iq_(-1),
        iu_(-1),
        m_(0),
        mx_(0),
        my_(0),
        mz_(0),
        gx_(0),
        gy_(0),
        gz_(0),
        coarse_level_(0),
        ir_loop_3_(-1),
        ir_loop_9_(-1){};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverBiCgStab);

  /// Charm++ PUP::able migration constructor
  EnzoSolverBiCgStab(CkMigrateMessage* m)
      : Solver(m),
        A_(NULL),
        is_alpha_(-1),
        is_beta_n_(-1),
        is_beta_d_(-1),
        is_rho0_(-1),
        is_err_(-1),
        is_err0_(-1),
        is_err_min_(-1),
        is_err_max_(-1),
        is_omega_(-1),
        is_omega_n_(-1),
        is_omega_d_(-1),
        is_rr_(-1),
        is_r0s_(-1),
        is_c_(-1),
        is_bs_(-1),
        is_xs_(-1),
        is_bnorm_(-1),
        is_vr0_(-1),
        is_ys_(-1),
        is_vs_(-1),
        is_us_(-1),
        is_qs_(-1),
        is_dot_sync_(-1),
        is_iter_(-1),
        res_tol_(0.0),
        index_precon_(-1),
        iter_max_(0),
        ir_(-1),
        ir0_(-1),
        ip_(-1),
        iy_(-1),
        iv_(-1),
        iq_(-1),
        iu_(-1),
        m_(0),
        mx_(0),
        my_(0),
        mz_(0),
        gx_(0),
        gy_(0),
        gz_(0),
        coarse_level_(0),
        ir_loop_3_(-1),
        ir_loop_9_(-1)

  {}

  /// Charm++ Pack / Unpack function
  void pup(PUP::er& p);

  /// Main solver entry routine
  virtual void apply(std::shared_ptr<Matrix> A, Block* block) throw();

  /// Type of this solver
  virtual std::string type() const { return "bicgstab"; }

  /// Projects RHS and sets initial vectors R, R0, and P
  void start_2(EnzoBlock* enzo_block, CkReductionMsg* msg) throw();

  /// Entry into BiCgStab iteration loop, begins refresh on P
  void loop_0a(EnzoBlock* enzo_block, CkReductionMsg*) throw();
  void loop_0b(EnzoBlock* enzo_block, CkReductionMsg*) throw();
  void loop_0(EnzoBlock* enzo_block) throw();

  /// First preconditioner solve
  void loop_2(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  void loop_25(EnzoBlock* enzo_block) throw();

  /// First matrix-vector product, begins DOT(V,R0) and projection of
  /// Y and V
  void loop_4(EnzoBlock* enzo_block) throw();

  /// Shifts Y and V, begins, first vector updates, begins refresh on Q
  void loop_6(EnzoBlock* enzo_block, CkReductionMsg*) throw();

  /// Second preconditioner solve, begins refresh on Y
  void loop_8(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  void loop_85(EnzoBlock* enzo_block) throw();

  /// Second matrix-vector product, begins DOT(U,U), DOT(U,Q) and
  /// projection of Y and U
  void loop_10(EnzoBlock* enzo_block) throw();

  /// Shifts Y and U, second vector updates, begins DOT(R,R) and
  /// DOT(R,R0)
  void loop_12(EnzoBlock* enzo_block, CkReductionMsg*) throw();

  /// Updates search direction, begins update on iteration counter
  void loop_14(EnzoBlock* enzo_block, CkReductionMsg*) throw();

  /// End the solve
  void end(EnzoBlock* enzo_block, int retval) throw();

  void dot_recv_parent(EnzoBlock*, int, long double*,
                       const std::vector<int>& is_array, int i_function,
                       int iter);
  void dot_recv_children(EnzoBlock*, int, long double*,
                         const std::vector<int>& is_array, int i_function);

protected:  // methods
  /// internal routine to handle actual start to solver
  void compute_(EnzoBlock* enzo_block) throw();

  /// Allocate temporary Fields
  void allocate_temporary_(Block* block) {
    Field field = block->data()->field();
    field.allocate_temporary(ir_);
    field.allocate_temporary(ir0_);
    field.allocate_temporary(ip_);
    field.allocate_temporary(iy_);
    field.allocate_temporary(iv_);
    field.allocate_temporary(iq_);
    field.allocate_temporary(iu_);
  }

  /// Dellocate temporary Fields
  void deallocate_temporary_(Block* block) {
    Field field = block->data()->field();
    field.deallocate_temporary(ir_);
    field.deallocate_temporary(ir0_);
    field.deallocate_temporary(ip_);
    field.deallocate_temporary(iy_);
    field.deallocate_temporary(iv_);
    field.deallocate_temporary(iq_);
    field.deallocate_temporary(iu_);
  }

  // Inner product methods

  void inner_product_(EnzoBlock*, int, long double*,
                      const std::vector<int>& isa, CkCallback callback,
                      int i_function);
  void dot_compute_tree_(EnzoBlock*, int, long double*,
                         const std::vector<int>& is_array, int i_function,
                         int iter);
  void dot_send_parent_(EnzoBlock*, int, long double*,
                        const std::vector<int>& is_array, int i_function,
                        int iter);
  void dot_send_children_(EnzoBlock*, int, long double*,
                          const std::vector<int>& is_array, int i_function);
  void dot_save_(EnzoBlock*, int, long double*,
                 const std::vector<int>& is_array);
  void dot_load_(EnzoBlock*, int, long double*,
                 const std::vector<int>& is_array);
  void dot_clear_(EnzoBlock*, int, const std::vector<int>& is_array);
  void dot_clear_(EnzoBlock*, int, long double*);
  void dot_done_(EnzoBlock*, int i_function, const char* file, int line);
  void dot_increment_(EnzoBlock*, int, const std::vector<int>& is_array,
                      long double* dot_block);

protected:
  inline long double& scalar_(Block* block, int i_scalar) {
    return *block->data()->scalar_long_double().value(i_scalar);
  }

  bool is_singular_() {
    return (A_->is_singular() && solve_type_ != solve_tree);
  }

  Sync& s_dot_sync_(EnzoBlock* block) {
    return *block->data()->scalar_sync().value(is_dot_sync_);
  }

  int& s_iter_(EnzoBlock* block) {
    return *block->data()->scalar_int().value(is_iter_);
  }

  /// Register all refresh phases
  void new_register_refresh_();

protected:  // attributes
  // NOTE: change pup() function whenever attributes change

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Corresponding ScalarData id's for solve_type == solve_tree
  int is_alpha_;
  int is_beta_n_;
  int is_beta_d_;
  int is_rho0_;
  int is_err_;
  int is_err0_;
  int is_err_min_;
  int is_err_max_;
  int is_omega_;
  int is_omega_n_;
  int is_omega_d_;
  int is_rr_;
  int is_r0s_;
  int is_c_;
  int is_bs_;
  int is_xs_;
  int is_bnorm_;
  int is_vr0_;
  int is_ys_;
  int is_vs_;
  int is_us_;
  int is_qs_;
  int is_dot_sync_;
  int is_iter_;

  /// Convergence tolerance on the relative residual
  double res_tol_;

  /// Preconditioner (-1 if none)
  int index_precon_;

  /// Maximum number of allowed BiCgStab iterations
  int iter_max_;

  /// BiCgStab vector id's
  int ir_;
  int ir0_;
  int ip_;
  int iy_;
  int iv_;
  int iq_;
  int iu_;

  /// Block field attributes
  int m_;             /// product mx_*my_*mz_ for convenience
  int mx_, my_, mz_;  /// total block size
  int gx_, gy_, gz_;  /// ghost zones

  /// The level of the tree solve if solve_type == solve_tree
  int coarse_level_;

  /// Refresh id's
  int ir_loop_3_;
  int ir_loop_9_;
};

#endif /* ENZO_ENZO_SOLVER_BICGSTAB_HPP */
