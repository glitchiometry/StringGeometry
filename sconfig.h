
#ifndef SCONFIG_H
#define SCONFIG_H
#include "basics.h"
#include "system.h"
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

//#define gsl_minimizer_type_mm gsl_multimin_fdfminimizer_steepest_descent
//#define gsl_minimizer_type_mm gsl_multimin_fdfminimizer_conjugate_pr
#define gsl_minimizer_type_mm gsl_multimin_fdfminimizer_vector_bfgs2

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
void string_config_randomize_mobile_coords(string_config *sc, double epsilon);
void string_config_randomize_coords(string_config *sc, double epsilon);
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
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
} total_length_pars;

void total_length_pars_init(total_length_pars *tlp, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map);

// RESUME: change length_func_fdf and Hess_length_func to compute upper/lower entries only (and
//    symmetrize recursively) in dim x dim blocks. This should improve performance by a factor between 2 and 4
//    depending on the dimension.
char grad_length_df_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H);
char grad_length_fdf_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H, gsl_vector *incr_f);
int grad_length_fdf(const gsl_vector *x, void *pars, gsl_vector *df, gsl_matrix *H);
int grad_length_df(const gsl_vector *x, void *pars, gsl_matrix *H);
int grad_length_f(const gsl_vector *x, void *pars, gsl_vector *df);

double mobile_length_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq);
double mobile_length(const gsl_vector *c_data, void *pars);

double total_length_f_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq);
double total_length_f(const gsl_vector *c_data, void *pars);
void total_length_df(const gsl_vector *x, void *pars, gsl_vector *df);
void total_length_df_exp(const gsl_vector *x, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, gsl_vector *df); // RESUME: consider passing tlf_pars instead of sc, c_map, top, etc.
void total_length_fdf_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, double *f, gsl_vector *df);
void total_length_fdf(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df);

void sc_solver_expand_cluster_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, int ci);
void sc_solver_expand(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data, gsl_vector *nc_data);

void read_mobile_coords_string_config(string_config *sc, contr_nbrlist *top, array_int *cmobile, gsl_vector *c_data, array_int *c_map);
void write_mobile_coords_string_config(string_config *sc, contr_nbrlist *top, array_int *cmobile, gsl_vector *c_data, array_int *c_map);

// LENGTH MINIMIZERS

void sc_minimizer_init_common_exp(contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, total_length_pars *tlf_pars, string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0);
void sc_minimizer_set_clusters_exp(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double epsilon);
void sc_minimizer_set_clusters_exp2(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, gsl_vector *f_data, double merge_ratio, double prec);
int sc_minimizer_test_stability_sd(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol);
int sc_minimizer_test_stability_bfgs2(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol);
void sc_minimizer_check_merging(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, double merge_rad);

typedef struct
{
  string_config *sc;
  gsl_multiroot_fsolver *solver;
  gsl_multiroot_function solver_data;
  int n_vars;
  void *pars;
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
  array_int *c_map;
  gsl_vector *c_data;
} sc_minimizer_rf;

void sc_minimizer_rf_init_exp(sc_minimizer_rf *scf, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *x, int (*f)(const gsl_vector *, void *, gsl_vector *), void *pars);
void free_sc_minimizer_rf(sc_minimizer_rf *scf);
int sc_minimizer_rf_relax(sc_minimizer_rf *scm, double tol);
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
  //  total_length_pars tlf_pars;
  void *pars;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  double char_length;
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
  double core_radsq;
  double *ext_f;
  int ext_f_site;
  double L0;
  double L0_cfxd;
  double offset;
} heuristic_pars;

void heuristic_pars_init(heuristic_pars *hpars, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, int ext_f_site, double *ext_f, double L0, double k);

double heuristic_f(const gsl_vector *c_data, void *pars);
void heuristic_df(const gsl_vector *c_data, void *pars, gsl_vector *df);
void heuristic_fdf(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df);

int grad_heuristic_f(const gsl_vector *c_data, void *pars, gsl_vector *f);
int grad_heuristic_df(const gsl_vector *c_data, void *pars, gsl_matrix *df);
int grad_heuristic_fdf(const gsl_vector *c_data, void *pars, gsl_vector *f, gsl_matrix *df);

void sc_minimizer_mm_init(sc_minimizer_mm *scm, string_config *sc);
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
  gsl_vector *f_data;
  contr_nbrlist top;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
  //  total_length_pars tlf_pars;
  int n_vars;
  double char_length;
  char set_flag;
} sc_minimizer_sd;

void sc_minimizer_sd_total_length_pars(sc_minimizer_sd *sd, total_length_pars *tlfpars);
void sc_minimizer_sd_heuristic_pars(sc_minimizer_sd *scm, heuristic_pars *hpars, int ext_f_site, double *ext_f, double L0, double k);
void sc_minimizer_sd_set(sc_minimizer_sd *scm, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars, const gsl_vector *x);
void sc_minimizer_sd_init(sc_minimizer_sd *scm, string_config *sc);
void sc_minimizer_sd_init_exp(sc_minimizer_sd *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x);
void free_sc_minimizer_sd(sc_minimizer_sd *scm);
double sc_minimizer_sd_iterate(sc_minimizer_sd *scm, double dt, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars, double prec);
void sc_minimizer_sd_check_merging(sc_minimizer_sd *scm, double merge_ratio, double prec);
void sc_minimizer_sd_expand_clusters(sc_minimizer_sd *scm, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars, double prec);
int sc_minimizer_sd_test_clusters(sc_minimizer_sd *scm, double *dt, int N_steps, double prec, double merge_ratio, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars);
int sc_minimizer_sd_relax2next_cluster(sc_minimizer_sd *scm, double *dt, double merge_radius, double prec, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars);
int sc_minimizer_sd_relax(sc_minimizer_sd *scm, double *dt, double merge_radius, double prec, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars);
void sc_minimizer_sd_solve(sc_minimizer_sd *scm, double *dt, double merge_radius, double prec, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars);
void sc_minimizer_sd_relax_diag(sc_minimizer_sd *scm, double dt, int N_steps, FILE *ofile, void (*vf)(const gsl_vector *, void *, gsl_vector *), void *pars);
int sc_minimizer_sd_relax_fin(sc_minimizer_sd *scm, total_length_pars *tlf_pars, double dt, double merge_ratio, double prec);
int sc_minimizer_sd_solve_fin(sc_minimizer_sd *scm, total_length_pars *tlf_pars, double dt, double merge_ratio, double prec);
int sc_minimizer_sd_relax_hfin(sc_minimizer_sd *scm, heuristic_pars *h_pars, double dt, double merge_ratio, double prec);
int sc_minimizer_sd_solve_hfin(sc_minimizer_sd *scm, heuristic_pars *h_pars, double dt, double merge_ratio, double prec);

// Root polisher for minimization (using Newton's method or the hybrid algorithm)
typedef struct
{
  string_config *sc;
  gsl_vector *c_data;
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
  array_int *c_map;
  total_length_pars tlf_pars;
  int n_vars;
  gsl_multiroot_fdfsolver *solver;
  gsl_multiroot_function_fdf solver_data;
} sc_minimizer_fin;

void sc_minimizer_fin_init_exp(sc_minimizer_fin *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x);
void free_sc_minimizer_fin(sc_minimizer_fin *scm);
int sc_minimizer_fin_iterate(sc_minimizer_fin *scm);
///void sc_minimizer_fin_check_merging(sc_minimizer_fin *scm, double rad);
int sc_minimizer_fin_relax(sc_minimizer_fin *scm, double prec);
void sc_minimizer_fin_relax_diag(sc_minimizer_fin *scm, int N_steps, FILE *ofile);

// Root polisher for minimization (using Newton's method or the hybrid algorithm)
typedef struct
{
  string_config *sc;
  gsl_vector *c_data;
  contr_nbrlist *top;
  array_int *cmobile;
  array_int *cmobile_map;
  array_int *c_map;
  heuristic_pars h_pars;
  int n_vars;
  gsl_multiroot_fdfsolver *solver;
  gsl_multiroot_function_fdf solver_data;
} sc_minimizer_hfin;

void sc_minimizer_hfin_init_exp(sc_minimizer_hfin *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, int ext_f_site, double *ext_f, double L0, double k);
void free_sc_minimizer_hfin(sc_minimizer_hfin *scm);
int sc_minimizer_hfin_iterate(sc_minimizer_hfin *scm);
int sc_minimizer_hfin_relax(sc_minimizer_hfin *scm, double prec);
void sc_minimizer_hfin_relax_diag(sc_minimizer_hfin *scm, int N_steps, FILE *ofile);

// Lagrange solver

typedef struct
{
  string_config *sc;
  char solver_type;
  int n_vars;
  int tau_addr;
  gsl_multiroot_fdfsolver *solver;
  gsl_multiroot_function_fdf solver_data;
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
} sc_Lagrange_solver;

// NOTE: eventually, this should probably be replaced/augmented with custom 
// 	routines involving sparse matrices (for string configurations where
// 	each vertex has O(1) neighbors.)
// NOTE: these functions should be tested
int Lagrange_solver_f(const gsl_vector *c_data, void *params, gsl_vector *f);
int Lagrange_solver_df(const gsl_vector *c_data, void *params, gsl_matrix *df);
int Lagrange_solver_fdf(const gsl_vector *c_data, void *params, gsl_vector *f, gsl_matrix *df);
char Lagrange_solver_df_compute_incr(gsl_matrix *incr, const double *xi, const double *xii, double *disp, int wt, double *lensq, double *delta_f_i, double tau);
void add_lower_triangular(gsl_matrix *dest, gsl_matrix *src);
void sub_lower_triangular(gsl_matrix *dest, gsl_matrix *src);

void sc_Lagrange_solver_set(sc_Lagrange_solver *scs, int ext_f_site, double *ext_f, double L, double tol);
void sc_Lagrange_solver_init(sc_Lagrange_solver *scs, string_config *sc, int ext_f_site, double *ext_f, double L);
void sc_Lagrange_solver_from_minimizer(sc_Lagrange_solver *scs, sc_minimizer_mm *scm, int ext_f_site, double *ext_f, double L);
void sc_Lagrange_solver_init_exp(sc_Lagrange_solver *scs, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *x, double tau, int ext_f_site, double *ext_f, double L);
void free_sc_Lagrange_solver(sc_Lagrange_solver *scs);

double Lagrange_solver_f_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data, gsl_vector *f);
double sc_solver_compute_total_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);
double sc_solver_compute_fixed_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);
double sc_solver_compute_mobile_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data);

int sc_Lagrange_solver_relax_diag(sc_Lagrange_solver *scs, int N_steps, FILE *ofile);
int sc_Lagrange_solver_relax(sc_Lagrange_solver *scs, double tol);

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
void check_consistency_ctop(string_config *sc, contr_nbrlist *top);
int get_sc_minimizer_mm_counter();
int get_sc_minimizer_rf_counter();
int get_sc_minimizer_sd_counter();
void get_sc_minimizer_counters(int *a, int *b);
int get_total_length_f_counter();
int get_total_length_df_counter();
int get_total_length_fdf_counter();
void reset_total_length_f_counter();
void reset_total_length_df_counter();
void reset_total_length_fdf_counter();
void reset_total_length_counters();

void reset_sc_minimizer_sd_counter();

// Remove a value from an 'embedding' of integers
void remove_embedding(array_int *emb, array_int *emb_inv, int emb_index);

// Contractions of string configs
void contract_leaves(string_config *sc, contr_nbrlist *aux, array_int *bdry, array_int *mobile, array_int *is_mobile);
void contract_elbows(string_config *sc, contr_nbrlist *aux, array_int *elbows, array_int *mobile, array_int *is_mobile);

// Prepare a 'contractable list' (or 'embedding')
void prep_contr_list(array_int *src, array_int *dest, array_int *dest_map, int map_size);

double sc_solver_distsq_exp(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_map1, array_int *c_map1, const gsl_vector *x1);
double euclid_distsq(const double *x0, const double *x1, int len);
double euclid_normsq(const double *x, int len);
double sc_minimizer_shortest_dist_mobile(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data);
double sc_minimizer_shortest_dist(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data);

typedef struct
{
  int dim;
  array_voidstar coords;
  array_char is_fxd;
  array_int mobile;
  array_int fxd;
  array_int fm_addr;
} tack_set;

void tack_set_init(tack_set *ts, int dim);
void add_fxd_pt_tack_set(tack_set *ts, double *x);
void add_mbl_pt_tack_set(tack_set *ts, double *x);
void free_tack_set(tack_set *ts, void (*free_coord_elem)(void *));
double tack_set_var_x(tack_set *ts);
double *tack_set_pos(tack_set *ts, int i);

typedef struct
{
  tack_set *ts;
  array_voidstar tops; // edge-weighted graphs
  array_voidstar embs; // embeddings
} linked_sc;

void linked_sc_init(linked_sc *lsc, tack_set *ts);
char linked_sc_connected(linked_sc *lsc);
void free_linked_sc(linked_sc *lsc);
void free_linked_sc_shallow(linked_sc *lsc);
void add2linked_sc(linked_sc *lsc, edge_wtd_graph *top, array_int *emb);

// RESUME: implement new (undefined) functions in sconfig.c, consider "refactoring"
//  to incorporate 'coord_pars' in other structures below.
typedef struct
{
  tack_set *ts;
  array_int map;
  aarray_int fibers;
  array_int cmobile;
  array_int cmobile_map;
  array_int c_map;
} coord_pars;

void coord_pars_from_tack_set(coord_pars *cpars, tack_set *ts);
void free_coord_pars(coord_pars *cpars);
void coord_pars_merge(coord_pars *cpars, int ci, int cj);
double contr_pointset_dist(const gsl_vector *x0, coord_pars *cpars0, const gsl_vector *x1, coord_pars *cpars1);
void transcribe_coord_pars(coord_pars *src, coord_pars *dest);
char coord_pars_corresp(coord_pars *a, coord_pars *b);

double lsc_string_length(void *ctop_, void *cemb_, coord_pars *cpars, const gsl_vector *x);

void lsc_solver_transcribe_tops(array_voidstar *ctops, array_voidstar *tops);
void lsc_solver_transcribe_embs(array_voidstar *cembs, aarray_int *cembs_map, array_voidstar *embs, int range_size);

typedef struct
{
  linked_sc *lsc;
  // Coordinate parameters (including a map from 'tack set' to clusters and its 'inverse')
  coord_pars *cpars;
  // Contracted topologies (as edge-weighted graphs) 
  array_voidstar *ctops;
  // Contracted embeddings (as array_ints)
  array_voidstar *cembs;
  // Inverse of embeddings (over contracted point-set, or the image of ts_map)
  aarray_int *cembs_map;
  double *Ls;
  // "Force" constant
  double k;
  // Vertex or string config
  int i;
  // If non-null, the force exerted at vertex 'i' (or its associated cluster)
  double *ext_f;
} lsc_heuristic_pars;

typedef struct
{
  linked_sc *lsc;
  int i;
  double *ext_f;
  gsl_vector *c_data;
  lsc_heuristic_pars hpars; // RESUME: define these types!
  gsl_multimin_fdfminimizer *hsolver;
  gsl_multimin_function_fdf hsolver_data;
  coord_pars cpars;
  // Contracted topologies (as edge-weighted graphs) 
  array_voidstar ctops;
  // Contracted embeddings (as array_ints)
  array_voidstar cembs;
  // Inverse of embeddings (over contracted point-set, or the image of ts_map)
  aarray_int cembs_map;
  double *Ls;
  double char_length;
} lsc_h_minimizer;

void lsc_h_minimizer_init(lsc_h_minimizer *lscm, linked_sc *lsc, int i, double *ext_f, double *Ls, double k);
void free_lsc_h_minimizer(lsc_h_minimizer *lscm);
int lsc_h_minimizer_iterate(lsc_h_minimizer *lscm);
void lsc_h_minimizer_check_merging(lsc_h_minimizer *lscm, double merge_rad);
void lsc_h_minimizer_expand_clusters(lsc_h_minimizer *lscm);
int lsc_h_minimizer_relax(lsc_h_minimizer *lscm, double tol);
int lsc_h_minimizer_solve(lsc_h_minimizer *lscm, double tol);

double lsc_h_f(const gsl_vector *x, void *pars);
void lsc_h_df(const gsl_vector *x, void *pars, gsl_vector *df);
void lsc_h_fdf(const gsl_vector *x, void *pars, double *f, gsl_vector *df);

// RESUME: implement these! (For polishing initial guesses, to avoid spurious topologies)
int grad_lsc_h_f(const gsl_vector *x, void *pars, gsl_vector *f);
int grad_lsc_h_df(const gsl_vector *x, void *pars, gsl_matrix *J);
int grad_lsc_h_fdf(const gsl_vector *x, void *pars, gsl_vector *f, gsl_matrix *J);

typedef struct
{
  linked_sc *lsc;
  coord_pars *cpars;
  // Contracted topologies (as edge-weighted graphs) 
  array_voidstar *ctops;
  // Contracted embeddings (as array_ints)
  array_voidstar *cembs;
  // Inverse of embeddings (over contracted point-set, or the image of ts_map)
  aarray_int *cembs_map;
  // Coordinate map for Lagrange multipliers
  array_int *l_map;
  double *Ls;
  int i;
  double *ext_f;
} lsc_Lagrange_pars;

void lsc_Lagrange_pars_init(lsc_Lagrange_pars *lpars, linked_sc *lsc, coord_pars *cpars, array_voidstar *ctops, array_voidstar *cembs, aarray_int *cembs_map, double *Ls, int i, double *ext_f, array_int *l_map);

typedef struct
{
  // Associated linked string config
  linked_sc *lsc;
  // An index (either a distinguished string config to be minimized, or a (mobile) vertex at which to apply the external force)
  int i;
  double *ext_f;
  lsc_Lagrange_pars lpars;
  gsl_multiroot_fdfsolver *Lsolver;
  gsl_multiroot_function_fdf Lsolver_data;
  // Consider replacing this with a generic 'coordinate data' structure
  gsl_vector *c_data;
  coord_pars *cpars;
  // Contracted topologies (as edge-weighted graphs) 
  array_voidstar *ctops;
  // Contracted embeddings (as array_ints)
  array_voidstar *cembs;
  // Inverse of embeddings (over contracted point-set, or the image of ts_map)
  aarray_int *cembs_map;
  double *Ls;
  // Coordinate map for Lagrange multipliers
  array_int l_map;
  double char_length;
} lsc_L_solver;

void lsc_L_solver_init(lsc_L_solver *lscs, lsc_h_minimizer *lsm);
void free_lsc_L_solver(lsc_L_solver *lscs);
int lsc_L_solver_iterate(lsc_L_solver *lscs);
int lsc_L_solver_relax(lsc_L_solver *lscs, double tol);
void lsc_L_solver_solve(lsc_L_solver *lscs, double tol);

void lsc_h_minimizer_from_lsc_L_solver(lsc_h_minimizer *lscm, lsc_L_solver *lscs, double k); // RESUME: implement this!

int lsc_L_f(const gsl_vector *x, void *pars, gsl_vector *f);
int lsc_L_df(const gsl_vector *x, void *pars, gsl_matrix *J);
int lsc_L_fdf(const gsl_vector *x, void *pars, gsl_vector *f, gsl_matrix *J);

#endif
