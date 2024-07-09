// RESUME: fix the 'robust' const_L_solver

// RESUME: define accelerator structure with precomputed sample points distributed 
// 	over the unit sphere for constant length surfaces (and over an indeterminate
// 	manifold in the minimum length case.)
#ifndef SCONFIG_H
#define SCONFIG_H
#include "basics.h"
#include "system.h"
#include "partition.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#define MIN_LEN_MODE 0
#define CONST_L_MODE 1
#define FMINIMIZER 0
#define FDFMINIMIZER 1
#define count_tally_size 256


typedef struct
{
  char mode;
  // Vertex coordinates (as pointers to externally allocated memory)
  array_voidstar pos;
  // External force (if any)
  double *ext_f;
  int ext_f_site;
  // Fixed/stationary vertices
  array_int fxd;
  array_char is_fxd;
  array_int fm_addr;
  array_int mobile;
  // Connectivity
  nbrlist top;
  // Edge weights
  aarray_int edge_wts;
  // Dimension
  int dim;
    // Number of spacial variables
  int n_s_vars;
} string_config;

void string_config_rotate(string_config *sc, gsl_matrix *Q, char mode);
void fprintf_string_config(string_config *sc, char *dirname);
void load_string_config(string_config *sc, array_voidstar *coords, char *dirname);
void string_line_coords(string_config *sc, double t, char mode);
int add_edge_string_config(string_config *sc, int i, int j, int wt);
void add_vertex_string_config(string_config *sc, double *x);

int string_config_init_loop(string_config *sc, array_int loop, int dim);
int string_config_init_seg(string_config *sc, array_int seq, int dim);
int string_config_init(string_config *sc, int nvs, int dim);
void free_string_config(string_config *sc);
void transcribe_string_config(string_config *src, string_config *dest);
void fix_point_string_config(string_config *sc, int i);
void unfix_point_string_config(string_config *sc, int i);
void set_pos_string_config(string_config *sc, int i, double *x);
void set_ext_force_string_config(string_config *sc, int i, double *f_ext_i);
void string_config_remove_ext_force(string_config *sc, int efi);
void string_config_line_coords(string_config *sc, double t);
// Consider preserving the original topology
void contract_string_config(string_config *sc, array_int *elem, array_int *cluster_map);
void expand_string_config(string_config *sc, int i, array_int *vs, nbrlist *local_top, aarray_int *local_edge_wts);
void string_config_compute_force(string_config *sc, double *sc_f);
void string_config_min_len(string_config *sc, double tol);
double string_config_min_len_step(string_config *sc, double tol, double stepsize);
int bisection_step_test(double *v, double *t, double *stepsize);

void string_config_centroid(string_config *sc, double *cntr);

double string_config_timestep(string_config *sc);
double string_config_sd_step_const_L(string_config *sc, double L, double tol, double stepsize); // RESUME: update this!
double string_config_total_length(string_config *sc);
double string_config_max_f(string_config *sc, double ref_val);
double *string_config_vertex_coords(string_config *sc, int i);
double string_config_var_x(string_config *sc);
const double *sc_solver_vertex_coords_exp(string_config *sc, int i, array_int *c_map, const gsl_vector *c_data);
const double *sc_solver_vertex_coords_exp2(string_config *sc, int ci, contr_nbrlist *top, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);
int string_config_relax_min_len(string_config *sc, double tol);
int string_config_relax_const_L(string_config *sc, int ext_f_site, double *ext_f, double L, double tol);
int string_config_relax_Lagrange(string_config *sc, int ext_f_site, double *ext_f, double L, double tol);

typedef struct
{
  string_config *sc;
  array_int *c_map;
  double core_radsq;
  partition *prt;
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
} total_length_func_pars;

// RESUME: change length_func_fdf and Hess_length_func to compute upper/lower entries only (and
//    symmetrize recursively) in dim x dim blocks. This should improve performance by a factor between 2 and 4
//    depending on the dimension.
void length_func_df_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H);
void length_func_fdf_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H, gsl_vector *incr_f);
int length_func_fdf(const gsl_vector *x, void *pars, gsl_vector *df, gsl_matrix *H);
int Hess_length_func(const gsl_vector *x, void *pars, gsl_matrix *H);
int grad_length_func(const gsl_vector *x, void *pars, gsl_vector *df);

double total_length_func_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq);
double total_length_func(const gsl_vector *c_data, void *pars);
void grad_total_length_func(const gsl_vector *x, void *pars, gsl_vector *df);
void grad_total_length_func_exp(const gsl_vector *x, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, gsl_vector *df); // RESUME: consider passing tlf_pars instead of sc, c_map, top, etc.
void comp_total_length_func_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, double *f, gsl_vector *df);
void comp_total_length_func(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df);

void sc_solver_expand_cluster_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, int ci);
void sc_solver_expand(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data, gsl_vector *nc_data);

void read_mobile_coords_string_config(string_config *sc, contr_nbrlist *top, array_int *cmobile, gsl_vector *c_data, array_int *c_map);
void write_mobile_coords_string_config(string_config *sc, contr_nbrlist *top, array_int *cmobile, gsl_vector *c_data, array_int *c_map);

// LENGTH MINIMIZERS

void sc_minimizer_init_common_exp(contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, total_length_func_pars *tlf_pars, string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0);
void sc_minimizer_set_clusters_exp(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double epsilon);
double sc_minimizer_var_mobile_coords(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *c_map, const gsl_vector *c_data);
int sc_minimizer_test_stability_sd(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol);
int sc_minimizer_test_stability_bfgs2(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol);
void sc_minimizer_check_merging(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, double merge_rad);

typedef struct
{
  string_config *sc;
  array_int c_map;
  gsl_multimin_fminimizer *solver;
  gsl_multimin_function solver_data;
  gsl_vector *smplx;
  double tol;
  gsl_vector *c_data;
  int n_vars;
  //  total_length_func_pars tlf_pars;
  void *pars;
  
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int cfxd;
  array_int cfxd_map;
  array_voidstar cfxd_coords;
} sc_minimizer_nm;

void sc_minimizer_nm_init(sc_minimizer_nm *scm, string_config *sc, double stepsize);
void free_sc_minimizer_nm(sc_minimizer_nm *scm);
int sc_minimizer_nm_iterate(sc_minimizer_nm *scm);
int sc_minimizer_nm_relax(sc_minimizer_nm *scm, double tol);
int sc_minimizer_nm_solve(sc_minimizer_nm *scm, double tol);
void sc_minimizer_nm_relax_diag(sc_minimizer_nm *scm, int n_steps, FILE *ofile);
//void sc_minimizer_nm_reset(sc_minimizer_nm *scf, double epsilon, const gsl_multimin_fminimizer_type *gsl_alg_type);
void sc_minimizer_nm_reset(sc_minimizer_nm *scf);
int get_sc_minimizer_nm_counter();
void sc_minimizer_nm_check_merging(sc_minimizer_nm *scm, double merge_rad);
void sc_minimizer_nm_init_const_L(sc_minimizer_nm *scm, string_config *sc, double L0, int ext_f_site, double *ext_f);

typedef struct
{
  string_config *sc;
  array_int c_map;
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function solver_data;
  double tol;
  gsl_vector *c_data;
  int n_vars;
  //  total_length_func_pars tlf_pars;
  void *pars;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
} sc_minimizer_rf;

int sc_minimizer_rf_relax(sc_minimizer_rf *scm, double tol);
void sc_minimizer_rf_init(sc_minimizer_rf *scf, string_config *sc, double prec);
void sc_minimizer_rf_reset(sc_minimizer_rf *scf, double epsilon, const gsl_multiroot_fsolver_type *gsl_alg_type); 
void free_sc_minimizer_rf(sc_minimizer_rf *scf);
int sc_minimizer_rf_iterate(sc_minimizer_rf *scf);
int sc_minimizer_rf_relax_diag(sc_minimizer_rf *scf, int n_steps, FILE *ofile);

typedef struct
{
  string_config *sc;
  array_int c_map;
  gsl_multimin_fdfminimizer *solver;
  gsl_multimin_function_fdf solver_data;
  gsl_vector *c_data;
  int n_vars;
  //  total_length_func_pars tlf_pars;
  void *pars;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int cfxd;
  array_int cfxd_map;
  array_voidstar cfxd_coords;
} sc_minimizer_mm;

// Heuristic function used to find initial guesses for the Lagrange solver
typedef struct
{
  // H(x) = 0.5 * k * (L(x) - L0)^2 - ext_f . (x_ef - x_ef_0)
  // dH(x) = k (L(x) - L0) dL(x) - ext_f . dx_ef
  string_config *sc;
  double k; // Steepness of constraining potential
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
  array_int *c_map;
  array_int *cfxd;
  array_int *cfxd_map;
  array_voidstar *cfxd_coords;
  double core_radsq;
  double *ext_f;
  int ext_f_site;
  double L0;
  //  double L0_fxd;
  double L0_cfxd;
} heuristic_pars;

double const_L_heuristic(const gsl_vector *c_data, void *pars);
void const_L_heuristic_df(const gsl_vector *c_data, void *pars, gsl_vector *df);
void const_L_heuristic_fdf(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df);

void sc_minimizer_mm_init(sc_minimizer_mm *scm, string_config *sc, double stepsize, double tol);
void free_sc_minimizer_mm(sc_minimizer_mm *scm);
int sc_minimizer_mm_iterate(sc_minimizer_mm *scm);
int sc_minimizer_mm_relax_diag(sc_minimizer_mm *scm, int n_steps, FILE *ofile);
int sc_minimizer_mm_relax(sc_minimizer_mm *scm, double tol);
int sc_minimizer_mm_solve(sc_minimizer_mm *scm, double tol);
void sc_minimizer_mm_contract(sc_minimizer_mm *scm, double epsilon, const gsl_multimin_fdfminimizer_type *gsl_alg_type, double stepsize, double tol); 
void sc_minimizer_mm_write_coords(sc_minimizer_mm *scm); // RESUME: implement this!
void sc_minimizer_mm_check_merging(sc_minimizer_mm *scm, double merge_rad);
void sc_minimizer_mm_reset(sc_minimizer_mm *scm);
void sc_minimizer_mm_init_exp(sc_minimizer_mm *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double stepsize, double tol);
void sc_minimizer_mm_init_heuristic(sc_minimizer_mm *scm, string_config *sc, int ext_f_site, double *ext_f, double L0, double k);
int sc_minimizer_mm_relax_heuristic(sc_minimizer_mm *scm, double tol);
int sc_minimizer_mm_solve_heuristic(sc_minimizer_mm *scm, double tol);
int sc_minimizer_mm_heuristic_update_L0_cmobile(sc_minimizer_mm *scm, heuristic_pars *hpars);

typedef struct
{
  string_config *sc;
  gsl_vector *c_data;
  gsl_vector *fdata;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
  //  total_length_func_pars tlf_pars;
  void *pars;
  int n_vars;
} sc_minimizer_sd;

void sc_minimizer_sd_init(sc_minimizer_sd *scm, string_config *sc);
void free_sc_minimizer_sd(sc_minimizer_sd *scm);
void sc_minimizer_sd_iterate(sc_minimizer_sd *scm, double dt);
void sc_minimizer_sd_check_merging(sc_minimizer_sd *scm, double rad);
int sc_minimizer_sd_relax(sc_minimizer_sd *scm, double dt, double merge_radius, double prec);
void sc_minimizer_sd_relax_diag(sc_minimizer_sd *scm, double dt, double merge_radius, int N_steps, FILE *ofile);
int get_sc_minimizer_sd_counter();

// Root polisher for minimization (using Newton's method or the hybrid algorithm)
typedef struct
{
  string_config *sc;
  gsl_vector *c_data;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
  array_int cfxd;
  array_int cfxd_map;
  array_voidstar cfxd_pos;
  total_length_func_pars tlf_pars;
  int n_vars;
  gsl_multiroot_fdfsolver *solver;
  gsl_multiroot_function_fdf solver_data;
} sc_minimizer_fin;

void sc_minimizer_fin_init(sc_minimizer_fin *scm, string_config *sc);
void sc_minimizer_fin_init_exp(sc_minimizer_fin *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x);
void free_sc_minimizer_fin(sc_minimizer_fin *scm);
int sc_minimizer_fin_iterate(sc_minimizer_fin *scm);
///void sc_minimizer_fin_check_merging(sc_minimizer_fin *scm, double rad);
int sc_minimizer_fin_relax(sc_minimizer_fin *scm, double prec);
void sc_minimizer_fin_relax_diag(sc_minimizer_fin *scm, int N_steps, FILE *ofile);

typedef struct
{
  string_config *sc;
  char solver_type;
  int n_vars;
  array_int c_map;
  gsl_multimin_fdfminimizer *solver;
  gsl_multimin_function_fdf solver_data;
  gsl_multimin_fminimizer *aux_solver;
  gsl_multimin_function aux_solver_data;
  double tol;
  gsl_vector *c_data;
  gsl_vector *smplx;
  //void *pars;
  total_length_func_pars tlf_pars;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  //double core_radsq; // RESUME: make sure this is initialized!
} sc_min_len_solver; // NOTE: this may have been rendered obsolete by sc_minimizer_mm and sc_minimizer_nm 

void sc_min_len_solver_cgmode(sc_min_len_solver *scs);
int sc_min_len_solver_n_vars(sc_min_len_solver *scs);
int sc_min_len_solver_relax(sc_min_len_solver *scs, double prec);
void sc_min_len_solver_relax_diag_simplex(sc_min_len_solver *scs, int N_steps, FILE *output);
void sc_min_len_solver_relax_diag_fdf(sc_min_len_solver *scs, int N_steps, FILE *output);
void sc_min_len_solver_init(sc_min_len_solver *scs, string_config *sc, double stepsize, double tol);
void free_sc_min_len_solver(sc_min_len_solver *scs);
void sc_min_len_solver_write_coords(sc_min_len_solver *scs);
void sc_min_len_solver_read_coords(sc_min_len_solver *scs, gsl_vector *c_data);
const double *sc_min_len_solver_vertex_coords(sc_min_len_solver *scs, int i, const gsl_vector *c_data);
gsl_vector *sc_min_len_solver_x(sc_min_len_solver *scs);


double const_L_target_func(const gsl_vector *c_data, void *pars);
//double const_L_target_func2(const gsl_vector *c_data, void *pars);
double const_L_target_func3(const gsl_vector *c_data, void *pars);

// RESUME: try to implement a variant of this that adjusts a (possibly random) Cartesian coordinate
//   to achieve constant length (rather than the forced vertex along the external force) and compare
//   performance for different choices (of coordinate/particle, as well as for various string configurations
//        and lengths.)

typedef struct
{
  string_config *sc;
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
  array_int *c_map;
  double L0;
  double *ext_f;
  int ext_f_site;
  double core_radsq;
} const_L_target_func_pars;

// Two step solver

typedef struct
{
  string_config *sc;
  double *ext_f;
  int ext_f_site;
  // Solver for minimizing total length over all other coordinates (treating ext_f_site as fixed)
  gsl_multimin_fdfminimizer *solver_a;
  gsl_multimin_function_fdf solver_a_data;
  // Solver for maximizing alignment with external force (treating all other vertices as fixed)
  gsl_multimin_fminimizer *solver_b;
  gsl_multimin_function solver_b_data;
  // Auxiliary data for solver_a (consider using an sc_minimizer_mm instance instead)
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
  
} sc_2step_solver;

// Lagrange solver

typedef struct
{
  string_config *sc;
  char solver_type;
  int n_vars;
  int tau_addr;
  gsl_multiroot_fdfsolver *solver;
  gsl_multiroot_function_fdf solver_data;
  gsl_multiroot_fsolver *fsolver;
  gsl_multiroot_function fsolver_data;
  double tol;
  gsl_vector *c_data;
  double mobile_length;
  double len_disc;
  double L;
  int ext_f_site;
  double *ext_f;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
  array_int cfxd;
  array_int cfxd_map;
  array_voidstar cfxd_coords;
} sc_Lagrange_solver;

// NOTE: eventually, this should probably be replaced/augmented with custom 
// 	routines involving sparse matrices (for string configurations where
// 	each vertex has O(1) neighbors.)
// NOTE: these functions should be tested
int Lagrange_solver_f(const gsl_vector *c_data, void *params, gsl_vector *f);
int Lagrange_solver_df(const gsl_vector *c_data, void *params, gsl_matrix *df);
int Lagrange_solver_fdf(const gsl_vector *c_data, void *params, gsl_vector *f, gsl_matrix *df);
void Lagrange_solver_df_compute_incr(gsl_matrix *incr, const double *xi, const double *xii, double *disp, int wt, double *lensq, double *delta_f_i, double tau);
void add_lower_triangular(gsl_matrix *dest, gsl_matrix *src);
void sub_lower_triangular(gsl_matrix *dest, gsl_matrix *src);

void sc_Lagrange_solver_set(sc_Lagrange_solver *scs, int ext_f_site, double *ext_f, double L, double tol); // RESUME: implement this! (The idea is to allow the total length and external force to be changed 'on the fly' for const_L solvers)
void sc_Lagrange_solver_init(sc_Lagrange_solver *scs, string_config *sc, int ext_f_site, double *ext_f, double L);
void sc_Lagrange_solver_from_minimizer(sc_Lagrange_solver *scs, sc_minimizer_mm *scm, int ext_f_site, double *ext_f, double L);
void free_sc_Lagrange_solver(sc_Lagrange_solver *scs);

double Lagrange_solver_f_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data, gsl_vector *f);
double sc_solver_compute_total_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);
double sc_solver_compute_fixed_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);
double sc_solver_compute_mobile_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);

int sc_Lagrange_solver_relax_f_diag(sc_Lagrange_solver *scs, int N_steps, FILE *ofile);
int sc_Lagrange_solver_relax_fdf_diag(sc_Lagrange_solver *scs, int N_steps, FILE *ofile);

int sc_Lagrange_solver_relax_fdf(sc_Lagrange_solver *scs, double tol);
int sc_Lagrange_solver_relax_f(sc_Lagrange_solver *scs, double tol);

// Composite solver used to iteratively find constant length configurations,
// refining a heuristic solver with each cycle.
typedef struct
{
  double tol;
  sc_Lagrange_solver scLs;
  sc_minimizer_mm scm;
  // An explicit pointer to scm.pars for simplified access
  heuristic_pars *hpars;
  // Reference coordinates for measuring convergence
  contr_nbrlist ref_top;
  array_int ref_cmobile;
  array_int ref_cmobile_map;
  array_int ref_c_map;
  gsl_vector *ref_data;
} sc_cl_composite_solver;

void sc_cl_composite_solver_init(sc_cl_composite_solver *sccs, string_config *sc, int ext_f_site, double *ext_f, double L, double tol);
void free_sc_cl_composite_solver(sc_cl_composite_solver *sccs);
double sc_cl_composite_solver_iterate(sc_cl_composite_solver *sccs);

// Utility functions
double rnd();
double min(double a, double b);
int fdf_iter_func(sc_min_len_solver *scs);
int f_iter_func(sc_min_len_solver *scs);

int get_sc_minimizer_mm_counter();
int get_sc_minimizer_rf_counter();
void get_sc_minimizer_counters(int *a, int *b);
int get_total_length_func_counter();
int get_grad_total_length_func_counter();
int get_comp_total_length_func_counter();

// Remove a value from an 'embedding' of integers
void remove_embedding(array_int *emb, array_int *emb_inv, int emb_index);

// Contractions of string configs
void group_leaves_elbows(string_config *sc, partition *prt);
void contract_leaves(contr_nbrlist *aux, array_int *bdry, array_int *mobile, array_int *is_mobile);
void contract_elbows(contr_nbrlist *aux, array_int *elbows, array_int *mobile, array_int *is_mobile);

// Prepare a 'contractable list' (or 'embedding')
void prep_contr_list(array_int *src, array_int *dest, array_int *dest_map, int map_size);

double sc_solver_distsq_exp(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_map1, array_int *c_map1, const gsl_vector *x1);
double euclid_distsq(const double *x0, const double *x1, int len);
double euclid_normsq(const double *x, int len);

double sc_minimizer_shortest_dist(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data);

#endif
