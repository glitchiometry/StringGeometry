#include "sconfig.h"
#define TRANSPOSE 1
#define NONTRANSPOSE 2
#define MAX_ITER 500
#define MAX_ITER_NM 1000
#define MAX_ITER_SD 10000

#define BRACKET_UPPER_BND 1e4
#define POWELL_SOLVER 2
#define BFGS2_SOLVER 0
#define CG_SOLVER 1
#define SIMPLEX 3

#define SC_CONST_L_SOLVER_ROTATED 1
#define SC_CONST_L_SOLVER_MINIMIZED 2

int total_length_func_counter = 0;
int grad_total_length_func_counter = 0;
int comp_total_length_func_counter = 0;

int sc_minimizer_mm_counter = 0;
int sc_minimizer_rf_counter = 0;
int sc_minimizer_nm_counter = 0;
int sc_minimizer_sd_counter = 0;

int get_sc_minimizer_sd_counter()
{
  return sc_minimizer_sd_counter;
}

int get_sc_minimizer_nm_counter()
{
  return sc_minimizer_nm_counter;
}

int get_sc_minimizer_mm_counter()
{
  return sc_minimizer_mm_counter;
}

int get_sc_minimizer_rf_counter()
{
  return sc_minimizer_rf_counter;
}

void get_sc_minimizer_counters(int *a, int *b)
{
  (*a) = sc_minimizer_mm_counter;
  (*b) = sc_minimizer_rf_counter;
}

int get_total_length_func_counter()
{
  return total_length_func_counter;
}

int get_grad_total_length_func_counter()
{
  return grad_total_length_func_counter;
}

int get_comp_total_length_func_counter()
{
  return comp_total_length_func_counter;
}

int glcs_count_tally[count_tally_size];

double matrix_get(gsl_matrix *A, int i, int j)
{
  return gsl_matrix_get(A, i, j);
}

double matrix_t_get(gsl_matrix *A, int i, int j)
{
  return gsl_matrix_get(A, j, i);
}

void line_coords(const double *x0, const double *v0, double t, double *xt, int len)
{
  for (int i = 0; i < len; i++) xt[i] = x0[i] + v0[i] * t;
}

double euclid_distsq(const double *x0, const double *x1, int len)
{
  double nsq = 0;
  for (int i = 0; i < len; i++)
    {
      double dxi = x1[i] - x0[i];
      nsq += dxi * dxi;
    }
  return nsq;
}

double euclid_normsq(const double *x, int len)
{
  double nsq = 0;
  for (int i = 0; i < len; i++) nsq += x[i] * x[i];
  return nsq;
}

void string_config_rotate(string_config *sc, gsl_matrix *Q, char mode)
{
  double (*get_op_Q)(gsl_matrix *, int, int) = mode != TRANSPOSE ? matrix_get : matrix_t_get;
  for (int i = 0; i < (*sc).pos.len; i++)
    {
      double *xi = (double *) (*sc).pos.e[i];
      double x0[(*sc).dim];
      for (int di = 0; di < (*sc).dim; di++) x0[di] = xi[di];
      for (int di = 0; di < (*sc).dim; di++)
	{
	  xi[di] = 0;
	  for (int dii = 0; dii < (*sc).dim; dii++) xi[di] += get_op_Q(Q, di, dii) * x0[dii];
	}
    }
  /*double ext_f[(*sc).dim];
  for (int di = 0; di < (*sc).dim; di++) ext_f[di] = (*sc).ext_f[di];
  for (int di = 0; di < (*sc).dim; di++)
    {
      (*sc).ext_f[di] = 0;
      for (int dii = 0; dii < (*sc).dim; dii++) (*sc).ext_f[di] += get_op_Q(Q, di, dii) * ext_f[dii];
      }*/
}

double string_config_timestep(string_config *sc)
{
  return 0.1 * sqrt(string_config_var_x(sc));
}

double *string_config_vertex_coords(string_config *sc, int i)
{
  return (double *) (*sc).pos.e[i];
}

void add_vertex_string_config(string_config *sc, double *x)
{
  int vi = (*sc).top.v.len;
  int init_nvars = (*sc).n_s_vars;
  add2array_voidstar(&((*sc).pos), x);
  extend_nbrlist(&((*sc).top));
  add2array_int(&((*sc).fm_addr), (*sc).mobile.len);
  add2array_int(&((*sc).mobile), vi);
  extend_aarray_int(&((*sc).edge_wts));
  (*sc).n_s_vars += (*sc).dim;
}

void remove_vertex_string_config(string_config *sc, int i)
{
  int fmi = (*sc).fm_addr.e[i];
  array_int *addrs;
  addrs = (*sc).is_fxd.e[i] ? &((*sc).fxd) : &((*sc).mobile);
  remove_array_int(addrs, fmi);
  if ((*addrs).len > 0)
    {
      int ii = (*addrs).e[fmi];
      (*sc).fm_addr.e[ii] = fmi;
    }
  remove_array_char(&((*sc).is_fxd), i);
  //if ((*sc).ext_f_site == i) (*sc).ext_f_site = -1;
  for (int ni = 0; ni < (*sc).top.v.e[i].len; ni++)
    {
      int ii = (*sc).top.v.e[i].e[ni];
      int i_i_ii = (*sc).top.i_of.e[i].e[ni];
      remove_array_int(&((*sc).edge_wts.e[ii]), i_i_ii);
    }
  remove_aarray_int(&((*sc).edge_wts), i);
  remove_vertex_nbrlist(&((*sc).top), i);
  remove_array_voidstar(&((*sc).pos), i, NULL);
}

int add_edge_string_config(string_config *sc, int i, int j, int wt)
{
  if (i != j)
    {
      // "'(i)ndex' of 'i' with respect to 'j'"
      int iji = add_edge_nbrlist_safe(&((*sc).top), i, j);
      if (iji > -1)
	{
	  int iij = (*sc).top.i_of.e[i].e[iji];
	  (*sc).edge_wts.e[i].e[iji] += wt;
	  (*sc).edge_wts.e[j].e[iij] += wt;
	}
      else
	{
	  add2array_int(&((*sc).edge_wts.e[i]), wt);
	  add2array_int(&((*sc).edge_wts.e[j]), wt);
	}
      return iji;
    }
  else
    {
      return -2;
    }
}

int string_config_init_loop(string_config *sc, array_int loop, int dim)
{
  printf("string_config_init_loop\n");
  int status = string_config_init_seg(sc, loop, dim);
  if (status == 0)
    {
      // Add the last edge from the loop
      add_edge_string_config(sc, loop.e[0], loop.e[loop.len - 1], 1);
    }
  printf("(done)\n");
  return status;
}

int string_config_init_seg(string_config *sc, array_int seq, int dim)
{
  printf("string_config_init_seg:\n");
  if (seq.len > 1 && dim > 1) {}
  else return seq.len < 2 | ((dim < 2) << 1);
  (*sc).dim = dim;
  int n_vertices = -1;
  for (int i = 0; i < seq.len; i++)
    {
      n_vertices = n_vertices >= seq.e[i] ? n_vertices : seq.e[i];
    }
  n_vertices += 1;
  string_config_init(sc, n_vertices, dim);
  int vim1 = seq.e[0];
  for (int i = 1; i < seq.len; i++)
    {
      int vi = seq.e[i];
      add_edge_string_config(sc, vim1, vi, 1);
      vim1 = vi;
    }
  printf("(done)\n");
  return 0;
}

// Initialize a string topology without coordinate assignments
int string_config_init(string_config *sc, int nvs, int dim)
{
  array_voidstar_init(&((*sc).pos), nvs);
  (*sc).pos.len = nvs;
  //(*sc).ext_f = (double *) calloc(dim, sizeof(double));
  array_int_init(&((*sc).fxd), 0);
  array_int_init(&((*sc).fm_addr), nvs);
  array_int_init(&((*sc).mobile), nvs);
  array_char_init(&((*sc).is_fxd), nvs);
  for (int i = 0; i < nvs; i++) 
    {
      (*sc).mobile.e[i] = i;
      (*sc).fm_addr.e[i] = i;
      (*sc).is_fxd.e[i] = 0;
    }
  (*sc).mobile.len = nvs;
  (*sc).fm_addr.len = nvs;
  nbrlist_init_precise(&((*sc).top), nvs);
  extend_nbrlist_n(&((*sc).top), nvs);
  aarray_int_init(&((*sc).edge_wts), nvs);
  (*sc).edge_wts.len = nvs;
  (*sc).dim = dim;
  return 0;
}

void free_string_config(string_config *sc)
{
  free_array_voidstar(&((*sc).pos), NULL);
  free_nbrlist(&((*sc).top));
  //  free((*sc).ext_f);
  free_array_int(&((*sc).fxd));
  free_array_int(&((*sc).mobile));
  free_array_int(&((*sc).fm_addr));
  free_array_char(&((*sc).is_fxd)); // RESUME: consider replacing this with a bit array
  free_aarray_int(&((*sc).edge_wts));
}

void transcribe_string_config(string_config *src, string_config *dest)
{
  transcribe_nbrlist(&((*src).top), &((*dest).top));
  transcribe_aarray_int(&((*src).edge_wts), &((*dest).edge_wts));
  transcribe_array_int(&((*src).fm_addr), &((*dest).fm_addr));
  transcribe_array_int(&((*src).fxd), &((*dest).fxd));
  transcribe_array_char(&((*src).is_fxd), &((*dest).is_fxd));
  transcribe_array_int(&((*src).mobile), &((*dest).mobile));
  array_voidstar_init(&((*dest).pos), (*src).pos.len);
  for (int i = 0; i < (*src).pos.len; i++)
    {
      double *xi = (double *) malloc(sizeof(double) * (*src).dim);
      double *xi0 = string_config_vertex_coords(src, i);
      for (int di = 0; di < (*src).dim; di++) xi[di] = xi0[di];
      add2array_voidstar(&((*dest).pos), xi0);
    }
}

void fix_point_string_config(string_config *sc, int i)
{
  if (!((*sc).is_fxd.e[i]))
    {
      int mi = (*sc).fm_addr.e[i];
      (*sc).fm_addr.e[i] = (*sc).fxd.len;
      (*sc).is_fxd.e[i] = 1;
      add2array_int(&((*sc).fxd), i);
      remove_array_int(&((*sc).mobile), mi); 
      if ((*sc).mobile.len > mi)
	{
	  int ii = (*sc).mobile.e[mi];
	  (*sc).fm_addr.e[ii] = mi;
	}
    }
}

void unfix_point_string_config(string_config *sc, int i)
{
  if ((*sc).is_fxd.e[i])
    {
      int addr_i = (*sc).fm_addr.e[i];
      remove_array_int(&((*sc).fxd), addr_i);
      if ((*sc).fxd.len > addr_i)
	{
	  int ii = (*sc).fxd.e[addr_i];
	  (*sc).fm_addr.e[ii] = addr_i;
	}
      int mi = (*sc).mobile.len;
      (*sc).fm_addr.e[i] = mi;
      add2array_int(&((*sc).mobile), i);
      (*sc).is_fxd.e[i] = 0;
    }
}

void set_pos_string_config(string_config *sc, int i, double *x)
{
  if (i < (*sc).pos.len) {}
  else
    {
      printf("Error (set_pos_string_config): attempting to set position of non-existent vertex %d of 0 - %d-1\n", i, (*sc).pos.len);
      exit(EXIT_FAILURE);
    }
  double *pos_i = string_config_vertex_coords(sc, i);
  if (pos_i != NULL)
    {
      for (int di = 0; di < (*sc).dim; di++) pos_i[di] = x[di];
    }
  else (*sc).pos.e[i] = x;
}

double string_config_var_x(string_config *sc)
{
  double vx = 0;
  double ax[(*sc).dim];
  for (int di = 0; di < (*sc).dim; di++) ax[di] = 0;
  int base_i = 0;
  int count = 0;
  for (int i = 0; i < (*sc).top.v.len; i++)
    {
      double *xi = string_config_vertex_coords(sc, i);
      if (xi != NULL) 
	{
	  count += 1;
	  for (int di = 0; di < (*sc).dim; di++) 
	    {
	      ax[di] += xi[di];
	      vx += xi[di] * xi[di];
	    }
	}
      base_i += (*sc).dim;
    }
  double inv_n_pts = 1. / count;
  double axsq = 0;
  for (int di = 0; di < (*sc).dim; di++) 
    {
      ax[di] *= inv_n_pts;
      axsq += ax[di] * ax[di];
    }
  return vx * inv_n_pts - axsq;
}

void string_config_compute_force_mobile_coords(string_config *sc, array_int *c_map, gsl_vector *x, gsl_vector *f)
{
  gsl_vector_set_zero(f);
  for (int mi = 0; mi < (*sc).mobile.len; mi++)
    {
      int i = (*sc).mobile.e[mi];
      double *xi = gsl_vector_ptr(x, (*c_map).e[mi]);
      double *fi = gsl_vector_ptr(f, (*c_map).e[mi]);
      for (int ni = 0; ni < (*sc).top.v.e[i].len; ni++)
	{
	  int ii = (*sc).top.v.e[i].e[ni];
	  if ((*sc).is_fxd.e[ii] || ii < i) {}
	  else continue;
	  double *xii;
	  if ((*sc).is_fxd.e[ii]) xii = string_config_vertex_coords(sc, ii);
	  else
	    {
	      int mii = (*sc).fm_addr.e[ii];
	      xii = gsl_vector_ptr(x, (*c_map).e[mii]);
	    }
	  double delx[(*sc).dim];
	  double delxsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] = xii[di] - xi[di];
	      delxsq += delx[di] * delx[di];
	    }
	  if (delxsq > 0) {}
	  else continue;
	  delxsq = (*sc).edge_wts.e[i].e[ni] / sqrt(delxsq);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] *= delxsq;
	      fi[di] += delx[di];
	    }
	  if (!(*sc).is_fxd.e[ii])
	    {
	      int mii = (*sc).fm_addr.e[ii];
	      double *fii = gsl_vector_ptr(f, (*c_map).e[mii]);
	      for (int di = 0; di < (*sc).dim; di++) fii[di] -= delx[di];
	    }
	}
    }
}

int string_config_max_total_weight(string_config *sc)
{
  int mtw = -1;
  for (int i = 0; i < (*sc).top.v.len; i++)
    {
      int lw = 0;
      for (int ni = 0; ni < (*sc).top.v.e[i].len; ni++)
	{
	  lw += (*sc).edge_wts.e[i].e[ni];
	}
      mtw = mtw >= lw ? mtw : lw;
    }
  return mtw;
}

void string_config_centroid(string_config *sc, double *cntr)
{
  for (int di = 0; di < (*sc).dim; di++) cntr[di] = 0;
  for (int i = 0; i < (*sc).top.v.len; i++)
    {
      double *xi = string_config_vertex_coords(sc, i);
      for (int di = 0; di < (*sc).dim; di++) cntr[di] += xi[di];
    }
  double wt = 1. / (*sc).top.v.len;
  for (int di = 0; di < (*sc).dim; di++) cntr[di] *= wt;
}

double string_config_total_length(string_config *sc)
{
  double L_ = 0;
  int n_vars = (*sc).n_s_vars;
  for (int i = 0; i < (*sc).top.v.len; i++)
    {
      double *pos_i = string_config_vertex_coords(sc, i);
      if (pos_i != NULL) {}
      else continue;
      for (int ni = 0; ni < (*sc).top.v.e[i].len; ni++)
	{
	  int ii = (*sc).top.v.e[i].e[ni];
	  if (ii < i) continue;
	  // Marginally faster would be to replace this with a binary operation, 
	  // and to store components in (partially empty) blocks of size 2^m
	  // (this would be slightly more memory efficient than using an aarray_double 
	  double *pos_ii = string_config_vertex_coords(sc, ii);
	  if (pos_ii != NULL) {}
	  else continue;
	  double distsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delxdi = pos_ii[di] - pos_i[di];
	      distsq += delxdi * delxdi;
	    }
	  L_ += (*sc).edge_wts.e[i].e[ni] * sqrt(distsq);
	}
    }
  return L_;
}

// RESUME: Check this!
// NOTE: it would probably be faster (and maybe even easier to check) just to
//        load each field independently. The main 'difficulty' would be choosing
//        between memory efficient files and 'debugging efficient' loading processes
//        when reading 'mobile', 'fxd', and 'is_fxd' arrays.
void load_string_config(string_config *sc, array_voidstar *coords, char *dirname)
{
  char fname[256];
  sprintf(fname, "%s/aux.dat", dirname);
  // Load 'metadata': the number of vertices, dimension, etc.
  int n_pts = -1;
  (*sc).dim = -1;
  FILE *ifile = fopen(fname, "r");
  if (ifile != NULL)
    {
      fscanf(ifile, "%d\n", &((*sc).dim));
      fscanf(ifile, "%d", &n_pts);
      fclose(ifile);
    }
  if (n_pts > -1)
    {
      //string_config_init(sc, n_pts, (*sc).dim);
      array_voidstar_init(coords, n_pts);
    }
  // Initialize the string configuration
  string_config_init(sc, n_pts, (*sc).dim);
  sprintf(fname, "%s/top.dat", dirname);
  char fname2[256];
  sprintf(fname2, "%s/wts.dat", dirname);
  ifile = fopen(fname, "r");
  FILE *ifile2 = fopen(fname2, "r");
  if (ifile != NULL && ifile2 != NULL)
    {
      int i = 0;
      while (1)
	{
	  int n_nbrs, n_wts;
	  int status = fscanf(ifile, "%d", &n_nbrs);
	  if (status == EOF) break;
	  fscanf(ifile2, "%d", &n_wts);
	  if (n_wts == n_nbrs) {}
	  else
	    {
	      printf("Something weird happened! oweiu928pu34442\n");
	      exit(EXIT_FAILURE);
	    }
	  for (int ni = 0; ni < n_nbrs; ni++)
	    {
	      int j, wt;
	      fscanf(ifile, "%d", &j);
	      fscanf(ifile2, "%d", &wt);
	      if (i > j) add_edge_string_config(sc, i, j, wt);
	    }
	  i += 1;
	}
      fclose(ifile);
      fclose(ifile2);
    }
  if ((*sc).dim > -1)
    {
      sprintf(fname, "%s/pos.dat", dirname);
      FILE *posfile = fopen(fname, "r");
      int i = 0;
      if (posfile != NULL)
	{
	  while (1)
	    {
	      double x0;
	      int status = fscanf(posfile, "%lg", &x0);
	      if (status == EOF) break;
	      double *xi = (double *) calloc((*sc).dim, sizeof(double));
	      xi[0] = x0;
	      for (int di = 1; di < (*sc).dim; di++)
		{
		  fscanf(posfile, "%lg", &(xi[di]));
		}
	      add2array_voidstar(coords, xi);
	      set_pos_string_config(sc, i, xi);
	      i += 1;
	    }
	  fclose(posfile);
	}
    }
  sprintf(fname, "%s/fxd.dat", dirname);
  ifile = fopen(fname, "r");
  if (ifile != NULL)
    {
      int i;
      while (fscanf(ifile, "%d", &i) != EOF)
	{
	  fix_point_string_config(sc, i);
	}
      fclose(ifile);
    }
}

// REQUIRES: system.h
void fprintf_string_config(string_config *sc, char *dirname)
{
  mkdir_s(dirname);
  char fname[256];
  sprintf(fname, "%s/aux.dat", dirname);
  FILE *ofile = fopen(fname, "w");
  if (ofile != NULL)
    {
      fprintf(ofile, "%d\n%d", (*sc).dim, (*sc).top.v.len);
      fclose(ofile);
    }
  sprintf(fname, "%s/top.dat", dirname);
  ofile = fopen(fname, "w");
  if (ofile != NULL)
    {
      fprintf_nbrlist(&((*sc).top), ofile);
      fclose(ofile);
    }
  sprintf(fname, "%s/pos.dat", dirname);
  ofile = fopen(fname, "w");
  if (ofile != NULL)
    {
      for (int i = 0; i < (*sc).top.v.len; i++)
	{
	  double *xi = string_config_vertex_coords(sc, i);
	  for (int di = 0; di < (*sc).dim; di++) fprintf(ofile, "%g ", xi[di]);
	  fprintf(ofile, "\n");
	}
      fclose(ofile);
    }
  sprintf(fname, "%s/wts.dat", dirname);
  ofile = fopen(fname, "w");
  if (ofile != NULL)
    {
      fprintf_aarray_int(&((*sc).edge_wts), ofile);
      fclose(ofile);
    }
  sprintf(fname, "%s/fxd.dat", dirname);
  ofile = fopen(fname, "w");
  if (ofile != NULL)
    {
      for (int fi = 0; fi < (*sc).fxd.len; fi++)
	{
	  fprintf(ofile, "%d ", (*sc).fxd.e[fi]);
	}
      fclose(ofile);
    }
}

double min(double a, double b)
{
  return a < b ? a : b;
}

void sc_solver_write_coords_exp(string_config *sc, const gsl_vector *x, array_int *c_map)
{
  for (int mi = 0; mi < (*sc).mobile.len; mi++)
    {
      const double *xi = gsl_vector_const_ptr(x, (*c_map).e[mi]);
      double *xi_ = string_config_vertex_coords(sc, (*sc).mobile.e[mi]);
      for (int di = 0; di < (*sc).dim; di++)
	{
	  xi_[di] = xi[di];
	}
    }
}

// RESUME: use this in the new minimizer_mm_init function
void sc_solver_read_mobile_coords_exp2(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_addr0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_addr1, array_int *c_map1, gsl_vector *x1)
{
  for (int cmi = 0; cmi < (*cmobile1).len; cmi++)
    {
      double *xi1 = gsl_vector_ptr(x1, (*c_map1).e[cmi]);
      int ci1 = (*cmobile1).e[cmi];
      int ri1 = (*top1).fibers.e[ci1].e[0];
      int ci0 = (*top0).map.e[ri1];
      int cmi0 = (*cmobile_addr0).e[ci0];
      const double *xi0;
      if (cmi0 > -1) xi0 = gsl_vector_const_ptr(x0, (*c_map0).e[cmi0]);
      else
	{
	  int ri0 = (*top0).fibers.e[ci0].e[0];
	  xi0 = string_config_vertex_coords(sc, ri0);
	}
      for (int di = 0; di < (*sc).dim; di++) xi1[di] = xi0[di];
    }
}

void sc_solver_read_coords_exp(string_config *sc, array_int *c_map, gsl_vector *c_data)
{
  for (int mi = 0; mi < (*sc).mobile.len; mi++)
    {
      int i = (*sc).mobile.e[mi];
      const double *xi = string_config_vertex_coords(sc, i);
      if (xi != NULL) {}
      else continue;
      double *xi_ = gsl_vector_ptr(c_data, (*c_map).e[mi]);
      for (int di = 0; di < (*sc).dim; di++) xi_[di] = xi[di];
    }
}

const double *sc_solver_vertex_coords_exp(string_config *sc, int i, array_int *c_map, const gsl_vector *c_data)
{
  return (*sc).is_fxd.e[i] ? string_config_vertex_coords(sc, i) : gsl_vector_const_ptr(c_data, (*c_map).e[(*sc).fm_addr.e[i]]);
}

const double *sc_solver_vertex_coords_exp2(string_config *sc, int ci, contr_nbrlist *top, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data)
{
  int i = (*top).fibers.e[ci].e[0];
  return (*cmobile_map).e[ci] == -1 ? string_config_vertex_coords(sc, i) : gsl_vector_const_ptr(c_data, (*c_map).e[(*cmobile_map).e[ci]]);
}

double mobile_length_func(const gsl_vector *c_data, void *pars)
{
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  return mobile_length_func_exp(c_data, (*tlf_pars).sc, (*tlf_pars).c_map, (*tlf_pars).top, (*tlf_pars).cmobile, (*tlf_pars).cmobile_map, (*tlf_pars).core_radsq);
}

double mobile_length_func_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq)
{
  double L = 0;
  nbrlist nw;
  aarray_int e_wts;
  array_int mobile;
  array_int mobile_map;
  if (top != NULL)
    {
      nw = (*top).top.top;
      e_wts = (*top).top.edge_wts;
      mobile = (*cmobile);
      mobile_map = (*cmobile_map);
    }
  else
    {
      //printf("Evaluating total_length_func with original topology and edge weights (exception!)\n");
      nw = (*sc).top;
      e_wts = (*sc).edge_wts;
      mobile = (*sc).mobile;
      mobile_map = (*sc).fm_addr;
    }
  for (int mci = 0; mci < mobile.len; mci++)
    {
      int ci = mobile.e[mci];
      const double *xi = gsl_vector_const_ptr(c_data, (*c_map).e[mci]);
      for (int ni = 0; ni < nw.v.e[ci].len; ni++)
	{
	  int cii = nw.v.e[ci].e[ni];
	  int rcii = cii;
	  
	  if (top != NULL)
	    {
	      rcii = (*top).fibers.e[cii].e[0];
	    }
	  if (cii < ci || (*sc).is_fxd.e[rcii]) {}
	  else continue;
	  const double *xii;
	  if ((*sc).is_fxd.e[rcii]) xii = string_config_vertex_coords(sc, rcii);
	  else xii = gsl_vector_const_ptr(c_data, (*c_map).e[mobile_map.e[cii]]);
	  double delxsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delx = xii[di] - xi[di];
	      delxsq += delx * delx;
	    }
	  if (delxsq <= core_radsq) continue;
	  double incr = e_wts.e[ci].e[ni] * sqrt(delxsq);
	  L += incr;
	}
    }
  //  printf("(done)\n");
  return L;  
}

double total_length_func_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq)
{
  //  printf("total_length_func_exp:\n");
  total_length_func_counter += 1;
  double L = 0;
  nbrlist nw;
  aarray_int e_wts;
  if (top != NULL)
    {
      nw = (*top).top.top;
      e_wts = (*top).top.edge_wts;
    }
  else
    {
      //printf("Evaluating total_length_func with original topology and edge weights (exception!)\n");
      nw = (*sc).top;
      e_wts = (*sc).edge_wts;
    }
  for (int ci = 0; ci < nw.v.len; ci++)
    {
      int ri;
      const double *xi;
      if (top != NULL)
	{
	  ri = (*top).fibers.e[ci].e[0];
	  xi = sc_solver_vertex_coords_exp2(sc, ci, top, cmobile_map, c_map, c_data);
	}
      else // RESUME: This should be phased out (solvers should probably use a contracted neighbor list)
	{
	  ri = ci;
	  xi = sc_solver_vertex_coords_exp(sc, ci, c_map, c_data);
	}
      for (int ni = 0; ni < nw.v.e[ci].len; ni++)
	{
	  int cii = nw.v.e[ci].e[ni];
	  if (cii < ci) {}
	  else continue;
	  int rii = cii;
	  const double *xii;
	  if (top != NULL)
	    {
	      rii = (*top).fibers.e[cii].e[0];
	      xii = sc_solver_vertex_coords_exp2(sc, cii, top, cmobile_map, c_map, c_data);
	    }
	  else
	    {
	      xii = sc_solver_vertex_coords_exp(sc, rii, c_map, c_data);
	    }
	  double delxsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delx = xii[di] - xi[di];
	      delxsq += delx * delx;
	    }
	  if (delxsq <= core_radsq) continue;
	  double incr = e_wts.e[ci].e[ni] * sqrt(delxsq);
	  L += incr;
	}
    }
  //  printf("(done)\n");
  return L;
}

double total_length_func(const gsl_vector *c_data, void *pars)
{
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  return total_length_func_exp(c_data, (*tlf_pars).sc, (*tlf_pars).c_map, (*tlf_pars).top, (*tlf_pars).cmobile, (*tlf_pars).cmobile_map, (*tlf_pars).core_radsq);
}

 void grad_total_length_func_exp(const gsl_vector *x, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, gsl_vector *df)
{
  grad_total_length_func_counter += 1;
  double *df_ = gsl_vector_ptr(df, 0);
  gsl_vector_set_zero(df);
  array_int mobile, mobile_addr;
  nbrlist nw;
  aarray_int e_wts;
  if (top != NULL && cmobile != NULL && cmobile_map != NULL)
    {
      nw = (*top).top.top;
      mobile = (*cmobile);
      mobile_addr = (*cmobile_map);
      e_wts = (*top).top.edge_wts;
    }
  else
    {
      nw = (*sc).top;
      mobile = (*sc).mobile;
      mobile_addr = (*sc).fm_addr;
      e_wts = (*sc).edge_wts;
    }
  for (int mi = 0; mi < mobile.len; mi++)
    {
      int ci = mobile.e[mi];
      int ri = ci;
      const double *xi;
      if (top != NULL)
	{
	  ri = (*top).fibers.e[ci].e[0];
	}
      xi = gsl_vector_const_ptr(x, (*c_map).e[mobile_addr.e[ci]]);
      double *fi = gsl_vector_ptr(df, (*c_map).e[mobile_addr.e[ci]]);
      for (int ni = 0; ni < nw.v.e[ci].len; ni++)
	{
	  int cii = nw.v.e[ci].e[ni];
	  int rii = cii;
	  if (top != NULL) rii = (*top).fibers.e[cii].e[0];
	  if ((rii < ri) || (*sc).is_fxd.e[rii]) {}
	  else continue;
	  double delx[(*sc).dim];
	  const double *xii;
	  double *fii;
	  if ((*sc).is_fxd.e[rii])
	    {
	      xii = string_config_vertex_coords(sc, rii);
	    }
	  else 
	    {
	      int mii = mobile_addr.e[cii];
	      xii = gsl_vector_const_ptr(x, (*c_map).e[mii]);
	      fii = gsl_vector_ptr(df, (*c_map).e[mii]);
	    }
	  if (xii != NULL) {}
	  else continue;
	  double delxsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] = xi[di] - xii[di];
	      delxsq += delx[di] * delx[di];
	    }
	  if (delxsq <= core_radsq)
	    {
	      continue;
	    }
	  delxsq = sqrt(delxsq);
	  delxsq = e_wts.e[ci].e[ni] / delxsq;
	  for (int di = 0; di < (*sc).dim; di++) 
	    {
	      delx[di] *= delxsq;
	      fi[di] += delx[di];
	    }
	  if (!(*sc).is_fxd.e[rii]) for (int di = 0; di < (*sc).dim; di++) fii[di] -= delx[di];
	}
    }
 }

/*
  L = sum_<i,j> w_ij |x_i-x_j|
  dL = sum_<i,j> w_ij (x_i-x_j)/|x_i-x_j|.(dx_i - dx_j)
  d^2L = sum_<i,j> w_ij ((dx_i-dx_j).(dx_i-dx_j)/|x_i-x_j| - (x_i-x_j).(dx_i-dx_j)(x_i-x_j).(dx_i-dx_j)/|x_i-x_j|^3)
  Elementary increments:
  w_ij delta_kl / |x_i-x_j|, H_ii, -H_ij, -H_ji, H_jj
  w_ij (x_i-x_j)(x_i-x_j) / |x_i-x_j|^3, -H_ii, -H_jj, H_ij, H_ji
  w_ij (x_i-x_j)/|x_i-x_j|, df_i, -df_j
  Note: DH_ii = w_ij (delta_kl - n_ij n_ij) / |x_i-x_j|, where df_i = w_ij n_ij
  Scaling: [sqrt]   vs.  [mult] * (dim - 1) dim / 2: direct approach is probably slightly faster in 2 and three dimensions.
  When evaluating both df and H:
          
 */

// NOTE: it could be slightly more efficient to compute increments as symmetric matrices (storing
//        upper or lower entries w/ diagonal) and then copying other entries for each mobile vertex
//    (and then symmetrizing the matrix). This is slightly more complicated, but has been implemented
//    successfully in the sc_Lagrange_solver.
void length_func_df_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H)
{
  int dim = (*incr_H).size1;
  double disp[dim];
  double delxsq = 0;
  for (int di = 0; di < dim; di++)
    {
      disp[di] = xi[di] - xj[di];
      delxsq += disp[di] * disp[di];
    }
  double inv_delxsq = 1. / delxsq;
  double inv_dist = sqrt(inv_delxsq);
  double factor = wt * inv_dist * inv_delxsq;
  gsl_matrix_set_zero(incr_H);
  for (int di = 0; di < dim; di++)
    {
      for (int dii = di; dii < dim; dii++)
	{
	  gsl_matrix_set(incr_H, di, dii, -disp[di] * disp[dii]);
	}
      gsl_matrix_set(incr_H, di, di, gsl_matrix_get(incr_H, di, di) + delxsq);
    }
  for (int di = 0; di < dim; di++)
    {
      for (int dii = 0; dii < di; dii++)
	{
	  gsl_matrix_set(incr_H, dii, di, gsl_matrix_get(incr_H, dii, di) * factor);
	  gsl_matrix_set(incr_H, di, dii, gsl_matrix_get(incr_H, dii, di));
	}
      gsl_matrix_set(incr_H, di, di, gsl_matrix_get(incr_H, di, di) * factor);
    }  
}

void length_func_fdf_incr(const double *xi, const double *xj, int wt, gsl_matrix *incr_H, gsl_vector *incr_f)
{
  gsl_vector_set_zero(incr_f);
  gsl_matrix_set_zero(incr_H);
  int dim = (*incr_H).size1;
  double disp[dim];
  double delxsq = 0;
  for (int di = 0; di < dim; di++)
    {
      disp[di] = xi[di] - xj[di];
      delxsq += disp[di] * disp[di];
    }
  double inv_dist = 1. / sqrt(delxsq);
  double factor = wt * inv_dist;
  for (int di = 0; di < dim; di++)
    {
      disp[di] *= inv_dist;
      gsl_vector_set(incr_f, di, disp[di] * wt);
    }
  gsl_matrix_set_identity(incr_H);
  for (int di = 0; di < dim; di++)
    {
      for (int dii = 0; dii < di; dii++)
	{
	  gsl_matrix_set(incr_H, di, dii, -disp[di] * disp[dii]);
	}
      gsl_matrix_set(incr_H, di, di, gsl_matrix_get(incr_H, di, di) - disp[di] * disp[di]);
    }
  for (int di = 0; di < dim; di++)
    {
      for (int dii = 0; dii < di; dii++)
	{
	  gsl_matrix_set(incr_H, di, dii, gsl_matrix_get(incr_H, di, dii) * factor);
	  gsl_matrix_set(incr_H, dii, di, gsl_matrix_get(incr_H, di, dii));
	}
      gsl_matrix_set(incr_H, di, di, gsl_matrix_get(incr_H, di, di) * factor);
    }
}

int length_func_fdf(const gsl_vector *x, void *pars, gsl_vector *df, gsl_matrix *H)
{
  //printf("length_func_fdf\n");
  if ((*df).size > 0) {}
  else return GSL_ENOPROGJ;
  gsl_vector_set_zero(df);
  gsl_matrix_set_zero(H);
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  contr_nbrlist *top = (*tlf_pars).top;
  string_config *sc = (*tlf_pars).sc;
  array_int *cmobile = (*tlf_pars).cmobile;
  array_int *cmobile_map = (*tlf_pars).cmobile_map;
  array_int *c_map = (*tlf_pars).c_map;
  double core_radsq = (*tlf_pars).core_radsq;
  //printf("test: %d\n", (*cmobile).len);
  // NOTE: it might be slightly faster to store increment matrices explicitly for each mobile vertex rather than have to allocate/deallocate these matrices with each function call (this approach would also be slightly easier to parallelize)
  gsl_matrix *incr = gsl_matrix_alloc((*sc).dim, (*sc).dim);
  gsl_vector *incr_f = gsl_vector_alloc((*sc).dim);
  for (int cmi = 0; cmi < (*cmobile).len; cmi++)
    {
      int ci = (*cmobile).e[cmi];
      const double *xci = gsl_vector_const_ptr(x, (*c_map).e[cmi]);
      gsl_vector_view dfci_ = gsl_vector_subvector(df, (*c_map).e[cmi], (*sc).dim);
      //double *dfci = gsl_vector_ptr(df, (*c_map).e[cmi]);
      gsl_matrix_view H_ci = gsl_matrix_submatrix(H, (*c_map).e[cmi], (*c_map).e[cmi], (*sc).dim, (*sc).dim);
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  //printf("%d ", cii);
	  if ((*cmobile_map).e[cii] < cmi) {}
	  else continue;
	  const double *xcii;
	  if ((*cmobile_map).e[cii] > -1) xcii = gsl_vector_const_ptr(x, (*c_map).e[(*cmobile_map).e[cii]]);
	  else
	    {
	      int fii = (*top).fibers.e[cii].e[0];
	      xcii = string_config_vertex_coords(sc, fii);
	    }
	  length_func_fdf_incr(xci, xcii, (*top).top.edge_wts.e[ci].e[ni], incr, incr_f);
	  gsl_matrix_add(&(H_ci.matrix), incr);
	  gsl_vector_add(&(dfci_.vector), incr_f);
	  if ((*cmobile_map).e[cii] > -1)
	    {
	      int cmii = (*cmobile_map).e[cii];
	      gsl_matrix_view H_cicii = gsl_matrix_submatrix(H, (*c_map).e[cmi], (*c_map).e[cmii], (*sc).dim, (*sc).dim);
	      gsl_matrix_view H_ciici = gsl_matrix_submatrix(H, (*c_map).e[cmii], (*c_map).e[cmi], (*sc).dim, (*sc).dim);
	      gsl_matrix_view H_cii = gsl_matrix_submatrix(H, (*c_map).e[cmii], (*c_map).e[cmii], (*sc).dim, (*sc).dim);
	      gsl_matrix_add(&(H_cii.matrix), incr);
	      gsl_matrix_sub(&(H_cicii.matrix), incr);
	      gsl_matrix_sub(&(H_ciici.matrix), incr);
	      gsl_vector_view f_cii = gsl_vector_subvector(df, (*c_map).e[cmii], (*sc).dim);
	      gsl_vector_sub(&(f_cii.vector), incr_f);
	    }
	}
      //printf("\n");
    }
  gsl_matrix_free(incr);
  gsl_vector_free(incr_f);
  //printf("(done)\n");
  return GSL_SUCCESS;
}

int Hess_length_func(const gsl_vector *x, void *pars, gsl_matrix *H)
{
  if ((*H).size1 > 0) {}
  else return GSL_ENOPROGJ;
  //printf("length_func_df:\n");
  gsl_matrix_set_zero(H);
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  contr_nbrlist *top = (*tlf_pars).top;
  string_config *sc = (*tlf_pars).sc;
  array_int *cmobile = (*tlf_pars).cmobile;
  array_int *cmobile_map = (*tlf_pars).cmobile_map;
  array_int *c_map = (*tlf_pars).c_map;
  double core_radsq = (*tlf_pars).core_radsq;
  // NOTE: it might be slightly faster to store increment matrices explicitly for each mobile vertex rather than have to allocate/deallocate these matrices with each function call (this approach would also be slightly easier to parallelize)
  gsl_matrix *incr = gsl_matrix_alloc((*sc).dim, (*sc).dim);
  for (int cmi = 0; cmi < (*cmobile).len; cmi++)
    {
      int ci = (*cmobile).e[cmi];
      const double *xci = gsl_vector_const_ptr(x, (*c_map).e[cmi]);
      //double *dfci = gsl_vector_ptr(df, (*c_map).e[cmi]);
      gsl_matrix_view H_ci = gsl_matrix_submatrix(H, (*c_map).e[cmi], (*c_map).e[cmi], (*sc).dim, (*sc).dim);
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  if ((*cmobile_map).e[cii] < cmi) {}
	  else continue;
	  const double *xcii;
	  if ((*cmobile_map).e[cii] > -1) xcii = gsl_vector_const_ptr(x, (*c_map).e[(*cmobile_map).e[cii]]);
	  else
	    {
	      int fii = (*top).fibers.e[cii].e[0];
	      xcii = string_config_vertex_coords(sc, fii);
	    }
	  length_func_df_incr(xci, xcii, (*top).top.edge_wts.e[ci].e[ni], incr);
	  gsl_matrix_add(&(H_ci.matrix), incr);
	  if ((*cmobile_map).e[cii] > -1)
	    {
	      int cmii = (*cmobile_map).e[cii];
	      gsl_matrix_view H_cicii = gsl_matrix_submatrix(H, (*c_map).e[cmi], (*c_map).e[cmii], (*sc).dim, (*sc).dim);
	      gsl_matrix_view H_ciici = gsl_matrix_submatrix(H, (*c_map).e[cmii], (*c_map).e[cmi], (*sc).dim, (*sc).dim);
	      gsl_matrix_view H_cii = gsl_matrix_submatrix(H, (*c_map).e[cmii], (*c_map).e[cmii], (*sc).dim, (*sc).dim);
	      gsl_matrix_add(&(H_cii.matrix), incr);
	      gsl_matrix_sub(&(H_cicii.matrix), incr);
	      gsl_matrix_sub(&(H_ciici.matrix), incr);
	    }
	}
    }
  gsl_matrix_free(incr);
  //printf("(done)\n");
  return GSL_SUCCESS;
}

int grad_length_func(const gsl_vector *x, void *pars, gsl_vector *df)
{
  if ((*df).size > 0) {}
  else return GSL_SUCCESS;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  grad_total_length_func_exp(x, (*tlf_pars).sc, (*tlf_pars).c_map, (*tlf_pars).top, (*tlf_pars).cmobile, (*tlf_pars).cmobile_map, (*tlf_pars).core_radsq, df);
  return GSL_SUCCESS;
}

void grad_total_length_func(const gsl_vector *x, void *pars, gsl_vector *df)
{
  if ((*df).size > 0) {}
  else return;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  grad_total_length_func_exp(x, (*tlf_pars).sc, (*tlf_pars).c_map, (*tlf_pars).top, (*tlf_pars).cmobile, (*tlf_pars).cmobile_map, (*tlf_pars).core_radsq, df);
}

 void comp_total_length_func_exp(const gsl_vector *c_data, string_config *sc, array_int *c_map, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, double core_radsq, double *f, gsl_vector *df)
{
  array_int mobile, mobile_addr;
  nbrlist nw;
  aarray_int e_wts;
  if (top != NULL && cmobile != NULL && cmobile_map != NULL)
    {
      mobile = (*cmobile);
      mobile_addr = (*cmobile_map);
      nw = (*top).top.top;
      e_wts = (*top).top.edge_wts;
    }
  else // RESUME: this case should be phased out
    {
      mobile = (*sc).mobile;
      mobile_addr = (*sc).fm_addr;
      nw = (*sc).top;
      e_wts = (*sc).edge_wts;
    }
  //printf("comp_total_length_func_exp:\n");
  comp_total_length_func_counter += 1;
  (*f) = 0;
  double *df_ = gsl_vector_ptr(df, 0);
  gsl_vector_set_zero(df);
  //  for (int i = 0; i < (*df).size; i++) df_[i] = 0;
  for (int ci = 0; ci < nw.v.len; ci++)
    {
      int ri = ci;
      const double *xi;
      if (top != NULL)
	{
	  ri = (*top).fibers.e[ci].e[0];
	  xi = sc_solver_vertex_coords_exp2(sc, ci, top, cmobile_map, c_map, c_data);
	}
      else xi = sc_solver_vertex_coords_exp(sc, ri, c_map, c_data);
      double *fi = NULL;
      if ((*sc).is_fxd.e[ri]) {}
      else
	{
	  int mi = mobile_addr.e[ci];
	  fi = &(df_[(*c_map).e[mi]]);
	}
      for (int ni = 0; ni < nw.v.e[ci].len; ni++)
	{
	  int cii = nw.v.e[ci].e[ni];
	  if (cii > ci) {}
	  else continue;
	  int rii = cii;
	  const double *xii;
	  if (top != NULL)
	    {
	      rii = (*top).fibers.e[cii].e[0];
	      xii = sc_solver_vertex_coords_exp2(sc, cii, top, cmobile_map, c_map, c_data);
	    }
	  else xii = sc_solver_vertex_coords_exp(sc, rii, c_map, c_data);
	  double *fii = NULL;
	  if ((*sc).is_fxd.e[rii]) {}
	  else
	    {
	      int mii = mobile_addr.e[cii];
	      fii = &(df_[(*c_map).e[mii]]);
	    }
	  double delx[(*sc).dim];
	  double delxsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] = xi[di] - xii[di];
	      delxsq += delx[di] * delx[di];
	    }
	  if (delxsq > core_radsq) {}
	  else continue;
	  double incr = e_wts.e[ci].e[ni] * sqrt(delxsq);
	  (*f) += incr;
	  if (fi != NULL || fii != NULL)
	    {
	      delxsq = incr / delxsq;
	      for (int di = 0; di < (*sc).dim; di++) delx[di] *= delxsq;
	      if (fi != NULL) for (int di = 0; di < (*sc).dim; di++) fi[di] += delx[di];
	      if (fii != NULL) for (int di = 0; di < (*sc).dim; di++) fii[di] -= delx[di];
	    }
	}
    }
}

void comp_total_length_func(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df)
{
  if ((*df).size > 0) {}
  else
  {
    (*f) = total_length_func(c_data, pars);
    return;
  }
  total_length_func_pars *tlf_pars = (total_length_func_pars *) pars;
  comp_total_length_func_exp(c_data, (*tlf_pars).sc, (*tlf_pars).c_map, (*tlf_pars).top, (*tlf_pars).cmobile, (*tlf_pars).cmobile_map, (*tlf_pars).core_radsq, f, df);  
} // END comp_total_length_func

// Heuristic function used to find initial guesses for the Lagrange solver
// RESUME: Determine what is wrong with this function (i.e. why it causes vertices to merge erroneously.)
double const_L_heuristic(const gsl_vector *c_data, void *pars)
{
  //  printf("const_L_heuristic:\n");
  heuristic_pars *hpars = (heuristic_pars *) pars;
  contr_nbrlist *top = (*hpars).top;
  array_int *cmobile = (*hpars).cmobile;
  array_int *cmobile_map = (*hpars).cmobile_map;
  array_int *c_map = (*hpars).c_map;
  double core_radsq = (*hpars).core_radsq;
  string_config *sc = (*hpars).sc;
  int c_efs = (*top).map.e[(*hpars).ext_f_site];
  int m_ef = (*cmobile_map).e[c_efs];
  const double *x_ef = gsl_vector_const_ptr(c_data, (*c_map).e[m_ef]);
  double H = 0;
  for (int di = 0; di < (*sc).dim; di++)
    {
      H -= x_ef[di] * (*hpars).ext_f[di];
    }
  double L_mobile = mobile_length_func_exp(c_data, sc, c_map, top, cmobile, cmobile_map, core_radsq); // RESUME: define this!
  double L = total_length_func_exp(c_data, sc, c_map, top, cmobile, cmobile_map, core_radsq);
  L -= (*hpars).L0;
  H += 0.5 * (*hpars).k * L * L;
  /*
  if (L > (*hpars).L0)
    {
      L -= (*hpars).L0;
      H += 0.5 * (*hpars).k * L * L;
    }
  */
  //  printf("(done: const_L_heuristic)\n");
  return H;
}

void const_L_heuristic_df(const gsl_vector *c_data, void *pars, gsl_vector *df)
{
  heuristic_pars *hpars = (heuristic_pars *) pars;
  contr_nbrlist *top = (*hpars).top;
  array_int *cmobile = (*hpars).cmobile;
  array_int *cmobile_map = (*hpars).cmobile_map;
  array_int *c_map = (*hpars).c_map;
  double core_radsq = (*hpars).core_radsq;
  string_config *sc = (*hpars).sc;
  int c_efs = (*top).map.e[(*hpars).ext_f_site];
  int m_ef = (*cmobile_map).e[c_efs];
  //const double *x_ef = gsl_vector_const_ptr(c_data, (*c_map).e[m_ef]);
  //  gsl_vector_set_zero(df); // RESUME: consider optimizing with a custom routine after verifying correctness  
  double L = total_length_func_exp(c_data, sc, c_map, top, cmobile, cmobile_map, core_radsq);
  //  if (L > (*hpars).L0)
  //{
  total_length_func_pars tlf_pars;
  tlf_pars.sc = sc;
  tlf_pars.top = top;
  tlf_pars.cmobile = cmobile;
  tlf_pars.cmobile_map = cmobile_map;
  tlf_pars.c_map = c_map;
  tlf_pars.core_radsq = core_radsq;
  comp_total_length_func(c_data, &tlf_pars, &L, df);
  double prefactor = (*hpars).k * (L - (*hpars).L0);
  gsl_vector_scale(df, prefactor);
  //}
  double *f_ef = gsl_vector_ptr(df, (*c_map).e[m_ef]);
  for (int di = 0; di < (*sc).dim; di++) f_ef[di] -= (*hpars).ext_f[di];
}

void const_L_heuristic_fdf(const gsl_vector *c_data, void *pars, double *f, gsl_vector *df)
{
  //printf("const_L_heuristic_fdf:\n");
  heuristic_pars *hpars = (heuristic_pars *) pars;
  contr_nbrlist *top = (*hpars).top;
  array_int *cmobile = (*hpars).cmobile;
  array_int *cmobile_map = (*hpars).cmobile_map;
  array_int *c_map = (*hpars).c_map;
  double core_radsq = (*hpars).core_radsq;
  string_config *sc = (*hpars).sc;
  int c_efs = (*top).map.e[(*hpars).ext_f_site];
  int m_ef = (*cmobile_map).e[c_efs];
  const double *x_ef = gsl_vector_const_ptr(c_data, (*c_map).e[m_ef]);
  //  gsl_vector_set_zero(df); // RESUME: consider optimizing with a custom routine after verifying correctness
  double L = total_length_func_exp(c_data, sc, c_map, top, cmobile, cmobile_map, core_radsq);
  //  if (L > (*hpars).L0)
  // {
  total_length_func_pars tlf_pars;
  tlf_pars.sc = sc;
  tlf_pars.top = top;
  tlf_pars.cmobile = cmobile;
  tlf_pars.cmobile_map = cmobile_map;
  tlf_pars.c_map = c_map;
  tlf_pars.core_radsq = core_radsq;
  grad_total_length_func(c_data, &tlf_pars, df);
  double diff = L - (*hpars).L0;
  double prefactor = (*hpars).k * diff;
  (*f) = prefactor * 0.5 * diff;
  gsl_vector_scale(df, prefactor);
  //}
  //else (*f) = 0;
  double *f_ef = gsl_vector_ptr(df, (*c_map).e[m_ef]);
  for (int di = 0; di < (*sc).dim; di++)
    {
      f_ef[di] -= (*hpars).ext_f[di];
      (*f) -= (*hpars).ext_f[di] * x_ef[di];
    }
  //printf("(done: const_L_heuristic_fdf)\n");
}

double *sc_minimizer_vertex_coords(string_config *sc, int i, contr_nbrlist *top, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data)
{
  int ci = (*top).map.e[i];
  int mci = (*cmobile_map).e[ci];
  if (mci > -1) return gsl_vector_ptr(c_data, (*c_map).e[mci]);
  else
    {
      int ri = (*top).fibers.e[ci].e[0];
      return string_config_vertex_coords(sc, ri);
    }
}

void sc_minimizer_mm_init_heuristic(sc_minimizer_mm *scm, string_config *sc, int ext_f_site, double *ext_f, double L0, double k)
{
  //  printf("sc_minimizer_mm_init_heuristic\n");
  (*scm).sc = sc;
  fix_point_string_config(sc, ext_f_site);
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  prep_contr_list(&((*sc).mobile), &((*scm).cmobile), &((*scm).cmobile_map), (*sc).top.v.len);
  //  printf("Test: %d %d\n", (*scm).top.top.top.v.len, (*scm).cmobile.len);
  // RESUME: fix memory error in contract_elbows
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  //printf("Contracting elbows\n");
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  unfix_point_string_config(sc, ext_f_site);
  int c_ef = (*scm).top.map.e[ext_f_site];
  (*scm).cmobile_map.e[c_ef] = (*scm).cmobile.len;
  add2array_int(&((*scm).cmobile), c_ef);
  //  printf("test 2: %d\n", (*scm).cmobile.len);
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  (*scm).c_map.len = (*scm).cmobile.len;
  heuristic_pars *hpars = (heuristic_pars *) malloc(sizeof(heuristic_pars));
  (*scm).pars = hpars;
  if ((*scm).cmobile.len > 0)
  {
    int n_vars = (*scm).cmobile.len * (*sc).dim;
    (*scm).n_vars = n_vars;
    (*scm).c_data = gsl_vector_alloc(n_vars);
    //(*scm).tol = stepsize;
    //printf("Setting solver_data\n");
    (*scm).solver_data.f = const_L_heuristic;
    (*scm).solver_data.df = const_L_heuristic_df;
    (*scm).solver_data.fdf = const_L_heuristic_fdf;
    (*scm).solver_data.params = (*scm).pars;
    (*scm).solver_data.n = n_vars;
    //printf("\tReading mobile coordinates\n");
    double stepsize = 0.01 * sqrt(string_config_var_x(sc));
    read_mobile_coords_string_config(sc, &((*scm).top), &((*scm).cmobile), (*scm).c_data, &((*scm).c_map));
    double *c_data_ = gsl_vector_ptr((*scm).c_data, 0);
    for (int di = 0; di < (*scm).n_vars; di++)
      {
	c_data_[di] += stepsize * (rnd() - 0.5);
      }
    for (int mci = 0; mci < (*scm).cmobile.len; mci++)
      {
	double *x_ci = gsl_vector_ptr((*scm).c_data, (*scm).c_map.e[mci]);
	for (int di = 0; di < (*sc).dim; di++) x_ci[di] += L0 * ext_f[di];
      }
    double *x_ef = sc_minimizer_vertex_coords(sc, ext_f_site, &((*scm).top), &((*scm).cmobile_map), &((*scm).c_map), (*scm).c_data);
    for (int di = 0; di < (*sc).dim; di++) x_ef[di] += L0 * ext_f[di];
    //printf("\t(done)\n");
    (*hpars).sc = sc;
    (*hpars).c_map = &((*scm).c_map);
    (*hpars).top = &((*scm).top);
    (*hpars).cmobile = &((*scm).cmobile);
    (*hpars).cmobile_map = &((*scm).cmobile_map);
    (*hpars).core_radsq = 1e-32;
    (*hpars).ext_f = ext_f;
    (*hpars).ext_f_site = ext_f_site;
    (*hpars).k = k;
    (*hpars).L0 = L0;
    //(*hpars).L0_fxd = string_config_fixed_length(sc);
    (*hpars).L0_cfxd = sc_solver_compute_fixed_length_exp(sc, &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), &((*scm).c_map), (*scm).c_data);
    (*hpars).offset = L0 - (*hpars).L0_cfxd;
    //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
    //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, n_vars);
    (*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
    //printf("Setting solver with stepsize %g\n", stepsize);
    gsl_multimin_fdfminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, stepsize, 1e-2);
    //printf("(done)\n");
  }
  else
  {
    printf("Something weird happened! a;lsdjfaos;ijdfaaw\n");
    // This case should not appear
    int n_vars = 0;
    (*scm).n_vars = n_vars;
    (*scm).c_data = NULL;
    (*hpars).sc = sc;
    (*hpars).c_map = &((*scm).c_map);
    (*hpars).top = &((*scm).top);
    (*hpars).cmobile = &((*scm).cmobile);
    (*hpars).cmobile_map = &((*scm).cmobile_map);
    (*hpars).core_radsq = 1e-32;
    //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
    //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, n_vars);
    (*scm).solver = NULL;
  }
}

int sc_minimizer_mm_relax_heuristic(sc_minimizer_mm *scm, double tol)
{
  printf("sc_minimizer_mm_relax_heuristic: TBD\n");
  return -1;
}

int sc_minimizer_mm_solve_heuristic(sc_minimizer_mm *scm, double tol)
{
  printf("sc_minimizer_mm_solve_heuristic: TBD\n");
  return -1;
}

int sc_minimizer_mm_heuristic_update_L0_cmobile(sc_minimizer_mm *scm, heuristic_pars *hpars)
{
  printf("sc_minimizer_mm_heuristic_update_L0_cmobile: TBD\n");
  return -1;
}

//const double k = 1.0;

double rnd()
{
  return ((double) rand()) / RAND_MAX;
}

double const_L_target_func3(const gsl_vector *c_data, void *pars)
{
  const_L_target_func_pars *cltf_pars = (const_L_target_func_pars *) pars;

  double L = total_length_func_exp(c_data, (*cltf_pars).sc, (*cltf_pars).c_map, (*cltf_pars).top, (*cltf_pars).cmobile, (*cltf_pars).cmobile_map, 1e-32);
  if ((*cltf_pars).L0 > L) {}
  else return L * sqrt(euclid_normsq(&((*cltf_pars).ext_f[0]), (*(*cltf_pars).sc).dim)) * 10;
  double dp = 0;
  int mi;
  if ((*cltf_pars).top != NULL) mi = (*(*cltf_pars).cmobile_map).e[(*(*cltf_pars).top).map.e[(*cltf_pars).ext_f_site]];
  else mi = (*(*cltf_pars).sc).fm_addr.e[(*cltf_pars).ext_f_site];
  const double *x_ef = gsl_vector_const_ptr(c_data, (*(*cltf_pars).c_map).e[mi]);
  for (int di = 0; di < (*(*cltf_pars).sc).dim; di++)
    {
      dp += x_ef[di] * (*cltf_pars).ext_f[di];
    }
  //printf("(done)\n");
  dp = -dp; // Consider removing the 0.01 * L term
  return dp;
}

double sc_solver_distsq_exp(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_map1, array_int *c_map1, const gsl_vector *x1)
{
  double distsq = 0;
  for (int mi = 0; mi < (*sc).mobile.len; mi++)
    {
      int i = (*sc).mobile.e[mi];
      int ci0 = (*top0).map.e[i];
      int cmi0 = (*cmobile_map0).e[ci0];
      const double *xi0;
      if (cmi0 > -1) xi0 = gsl_vector_const_ptr(x0, (*c_map0).e[cmi0]);
      else
	{
	  int ri0 = (*top0).fibers.e[ci0].e[0];
	  xi0 = string_config_vertex_coords(sc, ri0);
	}
      int ci1 = (*top1).map.e[i];
      int cmi1 = (*cmobile_map1).e[ci1];
      const double *xi1;
      if (cmi1 > -1) xi1 = gsl_vector_const_ptr(x1, (*c_map1).e[cmi1]);
      else
	{
	  int ri1 = (*top1).fibers.e[ci1].e[0];
	  xi1 = string_config_vertex_coords(sc, ri1);
	}
      distsq += euclid_distsq(xi0, xi1, (*sc).dim);
    }
  return distsq;
}

void read_mobile_coords_string_config(string_config *sc, contr_nbrlist *top, array_int *cmobile, gsl_vector *c_data, array_int *c_map)
{
  if (top != NULL && cmobile != NULL)
    {      
      //printf("Reading coordinates into contracted graph (cmobile.len = %d)\n", (*cmobile).len);
      int base = 0;
      for (int mi = 0; mi < (*cmobile).len; mi++)
	{
	  int ci = (*cmobile).e[mi];
	  //printf("cluster: %d (%d)\n", ci, (*top).top.top.v.len);
	  int i = (*top).fibers.e[ci].e[0];
	  //printf("Reading coordinates of %d (%d) into [%d,%d]\n", i, ci, base, base + (*sc).dim);
	  double *cdi = gsl_vector_ptr(c_data, base);
	  double *xi = string_config_vertex_coords(sc, i);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      cdi[di] = xi[di];
	    }
	  (*c_map).e[mi] = base;
	  base += (*sc).dim;
	}
    }
  else
    {
      int base = 0;
      for (int mi = 0; mi < (*sc).mobile.len; mi++)
	{
	  int ci = (*sc).mobile.e[mi];
	  double *cdi = gsl_vector_ptr(c_data, base);
	  double *xi = string_config_vertex_coords(sc, ci);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      cdi[di] = xi[di];
	    }
	  (*c_map).e[mi] = base;
	  base += (*sc).dim;
	}
      return;
    }  
}

void total_length_func_pars_init(total_length_func_pars *tlf_pars, string_config *sc, array_int *c_map, double epsilon)
{
  (*tlf_pars).core_radsq = epsilon * epsilon / (*sc).dim;
  (*tlf_pars).sc = sc;
  (*tlf_pars).c_map = c_map;
  (*tlf_pars).top = NULL;
  (*tlf_pars).cmobile = NULL;
  (*tlf_pars).cmobile_map = NULL;
}

// The idea is to use the steepest descent solver to test the stability
//  of clusters (until I can implement the 'Poisson' appproach.)
void sc_minimizer_sd_init(sc_minimizer_sd *scm, string_config *sc)
{
  (*scm).sc = sc;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) malloc(sizeof(total_length_func_pars));
  (*scm).pars = tlf_pars;
  (*tlf_pars).sc = sc;
  (*tlf_pars).top = &((*scm).top);
  (*tlf_pars).cmobile = &((*scm).cmobile);
  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
  (*tlf_pars).c_map = &((*scm).c_map);
  (*tlf_pars).core_radsq = 1e-16;
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  prep_contr_list(&((*sc).mobile), &((*scm).cmobile), &((*scm).cmobile_map), (*sc).top.v.len);
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  (*scm).c_map.len = (*scm).cmobile.len;
  int n_vars = (*scm).cmobile.len * (*sc).dim;
  (*scm).n_vars = n_vars;
  (*scm).c_data = gsl_vector_alloc(n_vars);
  (*scm).fdata = gsl_vector_alloc(n_vars);
  read_mobile_coords_string_config(sc, &((*scm).top), &((*scm).cmobile), (*scm).c_data, &((*scm).c_map));  
}

void free_sc_minimizer_sd(sc_minimizer_sd *scm)
{
  free((*scm).pars);
  free_array_int(&((*scm).c_map));
  gsl_vector_free((*scm).c_data);
  gsl_vector_free((*scm).fdata);
  free_array_int(&((*scm).cmobile));
  free_array_int(&((*scm).cmobile_map));
  free_contr_nbrlist(&((*scm).top));
}

void sc_minimizer_sd_iterate(sc_minimizer_sd *scm, double dt)
{
  sc_minimizer_sd_counter += 1;
  grad_total_length_func((*scm).c_data, (*scm).pars, (*scm).fdata);
  double *x = gsl_vector_ptr((*scm).c_data, 0);
  double *f = gsl_vector_ptr((*scm).fdata, 0);
  for (int i = 0; i < (*scm).c_data->size; i++)
    {
      x[i] -= f[i] * dt;
    }
}

void sc_minimizer_sd_check_merging(sc_minimizer_sd *scm, double rad)
{
  // NOTE: in general, this may require re-initializing (or 'resetting') a gsl_solver instance as well
  string_config *sc = (*scm).sc;
  contr_nbrlist *top = &((*scm).top);
  array_int *cmobile_map = &((*scm).cmobile_map);
  array_int *cmobile = &((*scm).cmobile);
  array_int *c_map = &((*scm).c_map);
  gsl_vector *c_data = (*scm).c_data;
  gsl_vector *fdata = (*scm).fdata;
  // void sc_minimizer_set_clusters_exp(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double epsilon) // RESUME: fix this!
  int nc0 = (*top).top.top.v.len;
  sc_minimizer_set_clusters_exp(sc, top, cmobile, cmobile_map, c_map, c_data, rad);
  int nc1 = (*top).top.top.v.len;
  if (nc0 != nc1)
    {
      // Reallocate coordinate and force data
      int n_vars = (*cmobile).len * (*sc).dim;
      gsl_vector *nc_data = gsl_vector_alloc(n_vars);
      gsl_vector *nfdata = gsl_vector_alloc(n_vars);
      (*scm).n_vars = n_vars;
      int base = 0;
      for (int ci = 0; ci < (*cmobile).len; ci++)
	{
	  double *xi = gsl_vector_ptr(nc_data, base);
	  double *fi = gsl_vector_ptr(nfdata, base);
	  double *xi0 = gsl_vector_ptr(c_data, (*c_map).e[ci]);
	  double *fi0 = gsl_vector_ptr((*scm).fdata, (*c_map).e[ci]);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      xi[di] = xi0[di];
	      fi[di] = fi0[di];
	    }
	  base += (*sc).dim;
	}
      gsl_vector_free(c_data);
      gsl_vector_free(fdata);
      (*scm).c_data = nc_data;
      (*scm).fdata = nfdata;
    }
}

int sc_minimizer_sd_relax(sc_minimizer_sd *scm, double dt, double merge_radius, double prec)
{
  string_config *sc = (*scm).sc;
  double r = 0.5;
  int prec_lim = log(prec / dt) / log(r) + 1;
  int count = 0;
  double timer_ = 0;
  double timer_lim = 5 * dt;
  int timer_counter = 0;
  double L_ref = total_length_func((*scm).c_data, (*scm).pars);
  gsl_vector *ave_pos = gsl_vector_alloc((*scm).n_vars);
  gsl_vector_set_zero(ave_pos);
  double *ave_pos_ = gsl_vector_ptr(ave_pos, 0);
  while (count < 10000000)
    {
      sc_minimizer_sd_iterate(scm, dt);
      timer_counter += 1;
      for (int ci = 0; ci < (*scm).n_vars; ci++)
	{
	  ave_pos_[ci] += gsl_vector_get((*scm).c_data, ci) * dt;
	}
      timer_ += dt;
      if (timer_ > timer_lim)
	{
	  double inv_time = 1. / timer_;
	  printf("inv_time = %g: ave_pos = ", inv_time);
	  for (int ci = 0; ci < (*scm).n_vars; ci++)
	    {
	      ave_pos_[ci] *= inv_time;
	      printf("%g ", ave_pos_[ci]);
	    }
	  printf("\n");
	  gsl_vector_memcpy((*scm).c_data, ave_pos);
	  printf("Evaluating rate of progress (count = %d)\n", count);
	  
	  
	  gsl_vector_set_zero(ave_pos);
	  timer_ = 0;
	}
      count += 1;
    }
  gsl_vector_free(ave_pos);
  return 1;
}

void sc_minimizer_sd_relax_diag(sc_minimizer_sd *scm, double dt, double merge_radius, int N_steps, FILE *ofile)
{
  for (int i = 0; i < N_steps; i++)
    {
      sc_minimizer_sd_iterate(scm, dt);
      sc_minimizer_sd_check_merging(scm, merge_radius);
      double *x = gsl_vector_ptr((*scm).c_data, 0);
      fprintf(ofile, "%d ", (*scm).cmobile.len);
      for (int mi = 0; mi < (*(*scm).sc).mobile.len; mi++)
	{
	  int i = (*(*scm).sc).mobile.e[mi];
	  int ci  = (*scm).top.map.e[i];
	  int mci = (*scm).cmobile_map.e[ci];
	  double *x_i;
	  if (mci > -1) x_i = &(x[(*scm).c_map.e[mci]]);
	  else x_i = string_config_vertex_coords((*scm).sc, i);
	  for (int di = 0; di < (*scm).sc->dim; di++)
	    {
	      fprintf(ofile, "%g ", x_i[di]);
	    }
	  fprintf(ofile, "\t");
	}
      //for (int ci = 0; ci < (*scm).c_data->size; ci++) fprintf(ofile, "%g ", x[ci]);
      double length = total_length_func((*scm).c_data, (*scm).pars);
      fprintf(ofile, "%g\n", length);
    }
}

int sc_minimizer_nm_relax(sc_minimizer_nm *scm, double tol)
{
  int count = 0;
  double merge_rad = 5 * tol;
  while (count < MAX_ITER_NM)
    {
      count += 1;
      int status;
      for (int i = 0; i < 5; i++) status = sc_minimizer_nm_iterate(scm);
      // Check if vertices need to be merged
      int n_vs0 = (*scm).n_vars;
      if (status == GSL_SUCCESS) sc_minimizer_nm_check_merging(scm, merge_rad);
      else
	{
	  printf("Iteration status = %d. Checking for merging vertices with default merge radius (1e-4)\n", status);
	  sc_minimizer_nm_check_merging(scm, 1e-4); // RESUME: update this (use geometry of configuration to define the merge radius)
	}
      if (n_vs0 == (*scm).n_vars)
	{
	  if (status == GSL_SUCCESS)
	    {
	      double size = gsl_multimin_fminimizer_size((*scm).solver);
	      status = gsl_multimin_test_size(size, tol);
	      if (status == GSL_CONTINUE) {}
	      else return status;
	    }
	  else return status;
	}
    }
}

// "Utility" relax functions (to be used to refactor 'solve' routines)
int scm_mm_relax_util(void *scm, double tol)
{
  return sc_minimizer_mm_relax((sc_minimizer_mm *) scm, tol);
}

int scm_nm_relax_util(void *scm, double tol)
{
  return sc_minimizer_nm_relax((sc_minimizer_nm *) scm, tol);
}

int scm_rf_relax_util(void *scm, double tol)
{
  return sc_minimizer_rf_relax((sc_minimizer_rf *) scm, tol);
}

// General purpose (utility/back-end) solver for minimizers
int sc_minimizer_solve_util(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **scm_solver_cdata, void *scm, int (*scm_relax_util)(void *, double))
{
  // RESUME: implement this to refactor explicit solvers (first, change 'relax' functions to transcribe solver->x to c_data
  //         after converging to a 'current best guess'.)
}

int sc_minimizer_nm_solve(sc_minimizer_nm *scm, double tol)
{
  //  printf("sc_minimizer_nm_solve\n");
  int count = 0;
  gsl_vector *x0 = gsl_vector_alloc((*scm).n_vars);
  gsl_vector_memcpy(x0, (*scm).solver->x);
  contr_nbrlist t0;
  transcribe_contr_nbrlist(&((*scm).top), &t0); 
  array_int cm0;
  transcribe_array_int(&((*scm).cmobile), &cm0);
  array_int cmm0;
  transcribe_array_int(&((*scm).cmobile_map), &cmm0);
  array_int cmap0;
  transcribe_array_int(&((*scm).c_map), &cmap0);
  double tolsq = tol * tol;
  int rstatus = GSL_CONTINUE;
  while (count < MAX_ITER)
    {
      count += 1;
      int status = sc_minimizer_nm_relax(scm, tol);
      // RESUME: check that status is acceptable before progressing (i.e. that the 'relaxation' hasn't terminated in a failure state)
      //double sc_solver_distsq_exp(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_map1, array_int *c_map1, const gsl_vector *x1)
      char realloc_flag = 0;
      if (contr_nbrlist_fiber_corresp(&t0, &((*scm).top), &((*(*scm).sc).mobile))) {}
      else
	{
	  printf("Fibers of relaxed contracted graph differ from those of initial graph\n");
	  status = GSL_CONTINUE;
	  realloc_flag = 1;
	}
      if (realloc_flag)
	{
	  // NOTE: there's probably a slightly more efficient approach to this (e.g. 'aligning' contraction t0 with (*scm).top)
	  // but this isn't (or shouldn't be) the rate-limiting step.
	  //printf("Reallocating t0, cm0, cmm0, cmap0, etc.\n");
	  free_contr_nbrlist(&t0);
	  free_array_int(&cm0);
	  free_array_int(&cmm0);
	  free_array_int(&cmap0);
	  gsl_vector_free(x0);
	  transcribe_contr_nbrlist(&((*scm).top), &t0);
	  transcribe_array_int(&((*scm).cmobile), &cm0);
	  transcribe_array_int(&((*scm).cmobile_map), &cmm0);
	  transcribe_array_int(&((*scm).c_map), &cmap0);
	  x0 = gsl_vector_alloc((*scm).n_vars);
	  //printf("(done)\n");
	}
      else
	{
	  double dispsq = sc_solver_distsq_exp((*scm).sc, &t0, &cm0, &cmm0, &cmap0, x0, &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), &((*scm).c_map), (*scm).solver->x);
	  if (dispsq > tolsq)
	    {
	      status = GSL_CONTINUE;
	    }
	  else
	    {
	      rstatus = GSL_SUCCESS;
	      break;
	    }
	}
      gsl_vector_memcpy(x0, (*scm).c_data);
      // RESET SOLVER AT CURRENT POSITION
      sc_minimizer_nm_reset(scm);
      
      // Consider perturbing mobile coordinates by random vector.
    }
  // Free auxiliary structures
  //printf("Freeing t0, cm0, cmm0, cmap0, x0\n");
  free_contr_nbrlist(&t0);
  free_array_int(&cm0);
  free_array_int(&cmm0);
  free_array_int(&cmap0);
  gsl_vector_free(x0);
  //printf("(done: sc_minimizer_nm_solve)\n");
  return rstatus;
}

void sc_minimizer_nm_init_const_L(sc_minimizer_nm *scm, string_config *sc, double L0, int ext_f_site, double *ext_f)
{
  fix_point_string_config(sc, ext_f_site);
  //printf("sc_minimizer_nm_init_const_L:\n");
  (*scm).sc = sc;
  const_L_target_func_pars *cltf_pars = (const_L_target_func_pars *) calloc(1, sizeof(const_L_target_func_pars));
  (*scm).pars = cltf_pars;
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  transcribe_array_int(&((*sc).mobile), &((*scm).cmobile));
  array_int_init(&((*scm).cmobile_map), (*sc).top.v.len);
  (*scm).cmobile_map.len = (*sc).top.v.len;
  for (int i = 0; i < (*sc).top.v.len; i++) (*scm).cmobile_map.e[i] = -1;
  for (int mi = 0; mi < (*scm).cmobile.len; mi++) (*scm).cmobile_map.e[(*scm).cmobile.e[mi]] = mi;
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  unfix_point_string_config(sc, ext_f_site);
  int c_ef = (*scm).top.map.e[ext_f_site];
  (*scm).cmobile_map.e[c_ef] = (*scm).cmobile.len;
  add2array_int(&((*scm).cmobile), c_ef);
  int n_vars = (*scm).cmobile.len * (*sc).dim;
  printf("test: %d %g %g\n", n_vars, string_config_total_length(sc), L0);
  (*scm).n_vars = n_vars;
  (*scm).solver = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, n_vars);
  //printf("test\n");
  (*scm).solver_data.f = const_L_target_func3;
  (*scm).solver_data.params = (*scm).pars;
  (*scm).solver_data.n = n_vars;
  (*cltf_pars).sc = sc;
  (*cltf_pars).top = &((*scm).top);
  (*cltf_pars).cmobile = &((*scm).cmobile);
  (*cltf_pars).cmobile_map = &((*scm).cmobile_map);
  (*cltf_pars).c_map = &((*scm).c_map);
  (*cltf_pars).core_radsq = 1e-32;
  (*cltf_pars).L0 = L0;
  (*cltf_pars).ext_f = ext_f;
  (*cltf_pars).ext_f_site = ext_f_site;
  (*scm).c_data = gsl_vector_alloc(n_vars);
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  (*scm).c_map.len = (*scm).cmobile.len;
  read_mobile_coords_string_config(sc, &((*scm).top), &((*scm).cmobile), (*scm).c_data, &((*scm).c_map));
  gsl_vector *smplx = gsl_vector_alloc(n_vars);
  double stepsize = sqrt(string_config_var_x(sc)) * 0.01;
  gsl_vector_set_all(smplx, stepsize);
  //printf("Setting solver\n");
  gsl_multimin_fminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, smplx);
  (*scm).smplx = smplx;
  //printf("(done)\n");
  //gsl_vector_free(smplx);
}

void sc_minimizer_nm_init(sc_minimizer_nm *scm, string_config *sc, double stepsize)
{
  (*scm).sc = sc;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) malloc(sizeof(total_length_func_pars));
  (*scm).pars = tlf_pars;
  // Determine if the topology is contractable
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  transcribe_array_int(&((*sc).mobile), &((*scm).cmobile));
  array_int_init(&((*scm).cmobile_map), (*sc).top.v.len);
  (*scm).cmobile_map.len = (*sc).top.v.len;
  for (int i = 0; i < (*sc).top.v.len; i++)
    {
      (*scm).cmobile_map.e[i] = -1;
    }
  for (int mi = 0; mi < (*scm).cmobile.len; mi++)
    {
      (*scm).cmobile_map.e[(*scm).cmobile.e[mi]] = mi;
    }
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  int n_vars = (*scm).cmobile.len * (*sc).dim;
  (*scm).n_vars = n_vars;
  (*scm).solver = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, n_vars);
  (*scm).solver_data.f = total_length_func;
  (*scm).solver_data.params = (*scm).pars;
  (*scm).solver_data.n = n_vars;
  (*tlf_pars).sc = sc;
  (*tlf_pars).top = &((*scm).top);
  (*tlf_pars).cmobile = &((*scm).cmobile);
  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
  (*tlf_pars).c_map = &((*scm).c_map);
  (*tlf_pars).core_radsq = 1e-32;
  (*scm).c_data = gsl_vector_alloc(n_vars);
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  (*scm).c_map.len = (*scm).cmobile.len;
  read_mobile_coords_string_config(sc, &((*scm).top), &((*scm).cmobile), (*scm).c_data, &((*scm).c_map));
  gsl_vector *smplx = gsl_vector_alloc(n_vars);
  gsl_vector_set_all(smplx, stepsize);
  gsl_multimin_fminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, smplx);
  (*scm).smplx = smplx;
  //gsl_vector_free(smplx);
}

void free_sc_minimizer_nm(sc_minimizer_nm *scm)
{
  free((*scm).pars);
  gsl_multimin_fminimizer_free((*scm).solver);
  free_array_int(&((*scm).c_map));
  gsl_vector_free((*scm).c_data);
  free_contr_nbrlist(&((*scm).top));
  free_array_int(&((*scm).cmobile));
  free_array_int(&((*scm).cmobile_map));
  gsl_vector_free((*scm).smplx);
}

// RESUME: test these after testing the characteristics of each solver with singular input
// (whether merging vertices converge to the same point, 'hover' around a common center (with some characteristic radius),
// or move unpredictably.)
void sc_minimizer_nm_reset(sc_minimizer_nm *scf)
{
  printf("sc_minimizer_nm_reset:\n");
  string_config *sc = (*scf).sc;
  int n_vars = (*(*scf).sc).dim * (*(*scf).sc).mobile.len;
  if (n_vars > (*scf).n_vars)
    {
      printf("Allocating new data for expanded coordinates\n");
      //gsl_vector_memcpy((*scf).c_data, (*scf).solver->x);
      gsl_vector *nc_data = gsl_vector_alloc(n_vars);
      sc_solver_expand((*scf).sc, &((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map), &((*scf).c_map), (*scf).solver->x, nc_data);
      (*scf).n_vars = n_vars;
      gsl_multimin_fminimizer_free((*scf).solver);
      gsl_vector_free((*scf).c_data);
      (*scf).c_data = nc_data;
      (*scf).solver = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, (*scf).n_vars);
      gsl_vector_free((*scf).smplx);
      gsl_vector *smplx = gsl_vector_alloc((*scf).n_vars);
      (*scf).solver_data.n = (*scf).n_vars; // RESUME: check if this step appears in mm_reset
      (*scf).smplx = smplx;
      printf("(done)\n");
    }
  double stepsize = sqrt(sc_minimizer_var_mobile_coords(sc, &((*scf).top), &((*scf).cmobile), &((*scf).c_map), (*scf).c_data)) * 0.1;
  gsl_vector_set_all((*scf).smplx, stepsize);
  gsl_multimin_fminimizer_set((*scf).solver, &((*scf).solver_data), (*scf).c_data, (*scf).smplx);
  printf("(done)\n");
}

void sc_minimizer_nm_check_merging(sc_minimizer_nm *scm, double merge_rad)
{
  //void sc_minimizer_check_merging(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, double merge_rad)
  int init_n_vs = (*scm).cmobile.len;
  contr_nbrlist *top = &((*scm).top);
  array_int *cmobile = &((*scm).cmobile);
  array_int *cmobile_map = &((*scm).cmobile_map);
  array_int *c_map = &((*scm).c_map);
  //gsl_vector_memcpy((*scm).c_data, (*scm).solver->x);
  gsl_vector *c_data = gsl_vector_alloc((*(*scm).c_data).size);
  gsl_vector_memcpy(c_data, (*scm).solver->x);
  sc_minimizer_check_merging((*scm).sc, top, cmobile, cmobile_map, c_map, &c_data, merge_rad);
  int final_n_vs = (*cmobile).len;
  if (final_n_vs != init_n_vs)
    {
      //      printf("sc_minimizer_nm_check_merging: resetting solver\n");
      // Reset the solver
      //printf("Resetting nm solver\n");
      (*scm).n_vars = (*c_data).size;
      //printf("Freeing solver\n");
      //      printf("Freeing simplex\n");
      gsl_vector_free((*scm).smplx);
      gsl_vector_free((*scm).c_data);
      //      printf("Freeing solver\n");
      gsl_multimin_fminimizer_free((*scm).solver); // check that this doesn't also free c_data
      //printf("Reallocating solver\n");
      (*scm).solver = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, (*scm).n_vars);
      //printf("Setting solver_data.n\n");
      (*scm).solver_data.n = (*c_data).size;
      //printf("Setting starting simplex\n");
      // double sc_minimizer_var_mobile_coords(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *c_map, const gsl_vector *c_data)
      gsl_vector *smplx = gsl_vector_alloc((*scm).n_vars);
      double stepsize = sqrt(sc_minimizer_var_mobile_coords((*scm).sc, top, cmobile, c_map, c_data)) / (*cmobile).len * 0.01;
      gsl_vector_set_all(smplx, stepsize);
      //printf("Setting solver\n");
      //      printf("Setting solver\n");
      (*scm).c_data = c_data;
      gsl_multimin_fminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, smplx);
      //printf("(done)\n");
      
      (*scm).smplx = smplx;
      //printf("(done)\n");
    }
}

int sc_minimizer_nm_iterate(sc_minimizer_nm *scm)
{
  sc_minimizer_nm_counter += 1;
  return gsl_multimin_fminimizer_iterate((*scm).solver);
}

void sc_minimizer_nm_relax_diag(sc_minimizer_nm *scm, int n_steps, FILE *ofile)
{
  if (ofile != NULL)
    {
      for (int i = 0; i < n_steps; i++)
	{
	  int status = sc_minimizer_nm_iterate(scm);
	  gsl_vector *x = (*(*scm).solver).x;
	  fprintf(ofile, "%d ", status);
	  for (int vi = 0; vi < (*x).size; vi++) fprintf(ofile, "%g ", gsl_vector_get(x, vi));
	  double l = total_length_func(x, (*scm).pars);
	  fprintf(ofile, ": %g %g\n", gsl_multimin_fminimizer_size((*scm).solver), l);
	}
    }
}

void check_consistency_cmobile(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map)
{
  for (int mi = 0; mi < (*cmobile).len; mi++)
    {
      int ci = (*cmobile).e[mi];
      if (ci < (*aux).top.top.v.len) {}
      else
	{
	  printf("Error: cmobile includes out of bounds values (%d:%d)\n", ci, (*aux).top.top.v.len);
	  exit(EXIT_FAILURE);
	}
      if ((*cmobile_map).e[ci] == mi) {}
      else
	{
	  printf("Error: inconsistency between cmobile_map and cmobile (%d %d: %d)\n", (*cmobile_map).e[ci], mi, ci);
	  exit(EXIT_FAILURE);
	}
    }
  for (int ci = 0; ci < (*aux).top.top.v.len; ci++)
    {
      if ((*cmobile_map).e[ci] > -1)
	{
	  int ri = (*aux).fibers.e[ci].e[0];
	  if ((*sc).is_fxd.e[ri])
	    {
	      printf("Error: cluster '%d' registered as mobile but contains a fixed vertex (%d)\n", ci, ri);
	      exit(EXIT_FAILURE);
	    }
	  if ((*cmobile).e[(*cmobile_map).e[ci]] == ci) {}
	  else
	    {
	      printf("Error: inconsistency in cmobile\n");
	      exit(EXIT_FAILURE);
	    }
	}
      else
	{
	  int ri = (*aux).fibers.e[ci].e[0];
	  if ((*sc).is_fxd.e[ri]) {}
	  else
	    {
	      printf("Error: cluster '%d' registered as fixed but is represented by a mobile vertex (%d)\n", ci, ri);
	      exit(EXIT_FAILURE);
	    }
	}
    }
}

void check_consistency_edge_wts(contr_nbrlist *aux)
{
  if ((*aux).top.top.v.len == (*aux).top.edge_wts.len && (*aux).fibers.len == (*aux).top.top.v.len)
    {
      for (int ci = 0; ci < (*aux).top.top.v.len; ci++)
	{
	  if ((*aux).top.top.v.e[ci].len == (*aux).top.edge_wts.e[ci].len) {}
	  else
	    {
	      printf("Error: discrepancy between contracted topology nbrlist[%d] of length %d and edge weights of length %d\n", ci, (*aux).top.top.v.e[ci].len, (*aux).top.edge_wts.e[ci].len);
	      exit(EXIT_FAILURE);
	    }
	}
    }
  else
    {
      printf("Error: discrepancy between total lengths of contracted topology (length %d), edge weights (length %d), or fibers (length %d)\n", (*aux).top.top.v.len, (*aux).top.edge_wts.len, (*aux).fibers.len);
    }
}

void sc_minimizer_init_common_exp(contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, total_length_func_pars *tlf_pars, string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0)
{
  // RESUME: change this to transcribe contracted topology and coordinate data explicitly
  contr_nbrlist_init(top, &((*sc).top), &((*sc).edge_wts));
  prep_contr_list(&((*sc).mobile), cmobile, cmobile_map, (*sc).top.v.len);
  contract_leaves(top, NULL, cmobile, cmobile_map);
  contract_elbows(top, NULL, cmobile, cmobile_map);
  array_int_init(c_map, (*cmobile).len);
  (*c_map).len = (*cmobile).len;
  int n_vars = (*cmobile).len * (*sc).dim;
  (*c_data) = gsl_vector_alloc(n_vars);
  sc_solver_read_mobile_coords_exp2(sc, top0, cmobile0, cmobile_map0, c_map0, x0, top, cmobile, cmobile_map, c_map, (*c_data));
  (*tlf_pars).sc = sc;
  (*tlf_pars).c_map = c_map;
  (*tlf_pars).top = top;
  (*tlf_pars).cmobile = cmobile;
  (*tlf_pars).cmobile_map = cmobile_map;
  (*tlf_pars).core_radsq = 1e-32;
}

void sc_minimizer_init_common(contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, total_length_func_pars *tlf_pars, string_config *sc)
{
  contr_nbrlist_init(top, &((*sc).top), &((*sc).edge_wts));
  // Contract elbows and leaves
  prep_contr_list(&((*sc).mobile), cmobile, cmobile_map, (*sc).top.v.len);
  contract_leaves(top, NULL, cmobile, cmobile_map);
  check_consistency_cmobile(sc, top, cmobile, cmobile_map);
  contract_elbows(top, NULL, cmobile, cmobile_map);
  check_consistency_edge_wts(top);
  check_consistency_cmobile(sc, top, cmobile, cmobile_map);
  array_int_init(c_map, (*cmobile).len);
  (*c_map).len = (*cmobile).len;
  int n_vars = (*cmobile).len * (*sc).dim;
  (*c_data) = gsl_vector_alloc(n_vars);
  read_mobile_coords_string_config(sc, top, cmobile, (*c_data), c_map);
  total_length_func_pars_init(tlf_pars, sc, c_map, 1e-32);
  (*tlf_pars).top = top;
  (*tlf_pars).cmobile = cmobile;
  (*tlf_pars).cmobile_map = cmobile_map;
}

void sc_minimizer_fin_init(sc_minimizer_fin *scf, string_config *sc)
{
  printf("sc_minimizer_fin_init:\n");
  sc_minimizer_init_common(&((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map), &((*scf).c_map), &((*scf).c_data), &((*scf).tlf_pars), sc);
  int n_vars = ((*scf).c_data -> size);
  (*scf).sc = sc;
  (*scf).solver_data.f = grad_length_func;
  (*scf).solver_data.df = Hess_length_func;
  (*scf).solver_data.fdf = length_func_fdf;
  (*scf).solver_data.params = &((*scf).tlf_pars);
  (*scf).solver_data.n = n_vars;
  (*scf).n_vars = n_vars;
  (*scf).solver = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_gnewton, (*scf).n_vars);
  gsl_multiroot_fdfsolver_set((*scf).solver, &((*scf).solver_data), (*scf).c_data);
  printf("(done)\n");
}

void sc_minimizer_fin_init_exp(sc_minimizer_fin *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x)
{
  printf("sc_minimizer_fin_init_exp:\n");
  // RESUME: Transcribe topology, cmobile, etc., define multiroot_function_fdf, set the solver
  printf("transcribing contr_nbrlist\n");
  transcribe_contr_nbrlist(top, &((*scm).top));
  transcribe_array_int(cmobile, &((*scm).cmobile));
  printf("cmobile: %d\n", (*scm).cmobile.len);
  transcribe_array_int(cmobile_map, &((*scm).cmobile_map));
  printf("cmobile_map: %d\n", (*scm).cmobile_map.len);
  transcribe_array_int(c_map, &((*scm).c_map));
  (*scm).c_data = gsl_vector_alloc((*x).size);
  gsl_vector_memcpy((*scm).c_data, x);
  (*scm).solver_data.f = grad_length_func;
  (*scm).solver_data.fdf = length_func_fdf;
  (*scm).solver_data.df = Hess_length_func;
  (*scm).solver_data.params = &((*scm).tlf_pars);
  (*scm).tlf_pars.top = top;
  (*scm).tlf_pars.cmobile = cmobile;
  (*scm).tlf_pars.cmobile_map = cmobile_map;
  (*scm).tlf_pars.c_map = c_map;
  (*scm).tlf_pars.core_radsq = 1e-32;
  (*scm).tlf_pars.sc = sc;
  int n_vars = (*x).size;
  (*scm).solver_data.n = n_vars;
  (*scm).n_vars = n_vars;
  (*scm).solver = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_gnewton, (*scm).n_vars);
  printf("Attempting to set solver\n");
  gsl_multiroot_fdfsolver_set((*scm).solver, &((*scm).solver_data), (*scm).c_data);
  printf("(done)\n");
}

void free_sc_minimizer_fin(sc_minimizer_fin *scm)
{
  free_contr_nbrlist(&((*scm).top));
  free_array_int(&((*scm).cmobile));
  free_array_int(&((*scm).cmobile_map));
  free_array_int(&((*scm).c_map));
  gsl_vector_free((*scm).c_data);
  gsl_multiroot_fdfsolver_free((*scm).solver);
}

int sc_minimizer_fin_iterate(sc_minimizer_fin *scm)
{
  return gsl_multiroot_fdfsolver_iterate((*scm).solver);
}

///void sc_minimizer_fin_check_merging(sc_minimizer_fin *scm, double rad);
int sc_minimizer_fin_relax(sc_minimizer_fin *scm, double prec)
{
  int rstatus = GSL_CONTINUE;
  int count = 0;
  while (count < MAX_ITER)
    {
      count += 1;
      int status = sc_minimizer_fin_iterate(scm);
      gsl_vector *dx = gsl_multiroot_fdfsolver_dx((*scm).solver);
      gsl_vector *f = gsl_multiroot_fdfsolver_f((*scm).solver);
      double *dx_ = gsl_vector_ptr(dx, 0);
      double *f_ = gsl_vector_ptr(f, 0);
      double err = sqrt(euclid_normsq(dx_, (*dx).size) + euclid_normsq(f_, (*dx).size));
      if (err < prec)
	{
	  rstatus = GSL_SUCCESS;
	  break;
	}
    }
}

void sc_minimizer_fin_relax_diag(sc_minimizer_fin *scm, int N_steps, FILE *ofile)
{
  for (int i = 0; i < N_steps; i++)
    {
      int status = sc_minimizer_fin_iterate(scm);
      if (ofile != NULL)
	{
	  fprintf(ofile, "%d ", status);
	  for (int di = 0; di < (*scm).n_vars; di++) fprintf(ofile, "%g ", gsl_vector_get((*scm).solver->x, di));
	  gsl_vector *dx = gsl_multiroot_fdfsolver_dx((*scm).solver);
	  gsl_vector *f = gsl_multiroot_fdfsolver_f((*scm).solver);
	  double len = total_length_func((*scm).solver->x, &((*scm).tlf_pars));
	  double err = sqrt(euclid_normsq(gsl_vector_ptr(dx, 0), (*dx).size) + euclid_normsq(gsl_vector_ptr(f, 0), (*dx).size));
	  fprintf(ofile, "%g %g\n", len, err);
	}
    }
}

// Root-finder minimizer (which may also be used for polishing)

void sc_minimizer_rf_init(sc_minimizer_rf *scf, string_config *sc, double stepsize)
{
  printf("sc_minimizer_rf_init:\n");
  (*scf).sc = sc;
  contr_nbrlist_init(&((*scf).top), &((*sc).top), &((*sc).edge_wts));
  // Contract elbows and leaves
  prep_contr_list(&((*sc).mobile), &((*scf).cmobile), &((*scf).cmobile_map), (*sc).top.v.len);
  contract_leaves(&((*scf).top), NULL, &((*scf).cmobile), &((*scf).cmobile_map));
  check_consistency_cmobile(sc, &((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map));
  printf("test\n");
  contract_elbows(&((*scf).top), NULL, &((*scf).cmobile), &((*scf).cmobile_map));
  check_consistency_edge_wts(&((*scf).top));
  check_consistency_cmobile(sc, &((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map));
  array_int_init(&((*scf).c_map), (*scf).cmobile.len);
  (*scf).c_map.len = (*scf).cmobile.len;
  int n_vars = (*scf).cmobile.len * (*sc).dim;
  (*scf).c_data = gsl_vector_alloc(n_vars);
  printf("Attempting to read mobile coords\n");
  read_mobile_coords_string_config(sc, &((*scf).top), &((*scf).cmobile), (*scf).c_data, &((*scf).c_map));
  printf("\t(done)\n");
  (*scf).solver_data.f = grad_length_func;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) malloc(sizeof(total_length_func_pars));
  total_length_func_pars_init(tlf_pars, sc, &((*scf).c_map), 1e-32);
  (*scf).pars = tlf_pars;
  (*tlf_pars).top = &((*scf).top);
  (*tlf_pars).cmobile = &((*scf).cmobile);
  (*tlf_pars).cmobile_map = &((*scf).cmobile_map);
  //  (*scf).tlf_pars.core_radsq = 1e-32;
  (*scf).solver_data.params = (*scf).pars;
  (*scf).solver_data.n = n_vars;
  (*scf).n_vars = n_vars;
  (*scf).solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, (*scf).n_vars);
  printf("Setting solver\n");
  gsl_multiroot_fsolver_set((*scf).solver, &((*scf).solver_data), (*scf).c_data);
  printf("(done)\n");
}

void free_sc_minimizer_rf(sc_minimizer_rf *scf)
{
  printf("free_sc_minimizer_rf\n");
  gsl_vector_free((*scf).c_data);
  free_array_int(&((*scf).c_map));
  gsl_multiroot_fsolver_free((*scf).solver);
  free_contr_nbrlist(&((*scf).top));
  free_array_int(&((*scf).cmobile));
  free_array_int(&((*scf).cmobile_map));
  printf("(done)\n");
}

// RESUME: update usage to include new arguments (c_map, cmobile)
double sc_minimizer_var_mobile_coords(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *c_map, const gsl_vector *c_data)
{
  double var = 0;
  double com[(*sc).dim];
  for (int i = 0; i < (*sc).dim; i++) com[i] = 0;
  for (int fi = 0; fi < (*sc).fxd.len; fi++)
    {
      int i = (*sc).fxd.e[fi];
      const double *xi = string_config_vertex_coords(sc, i);
      for (int di = 0; di < (*sc).dim; di++) com[di] += xi[di];
    }
  for (int mi = 0; mi < (*cmobile).len; mi++)
    {
      int ci = (*cmobile).e[mi];
      const double *xi = gsl_vector_const_ptr(c_data, (*c_map).e[mi]);
      if (top == NULL) for (int di = 0; di < (*sc).dim; di++) com[di] += xi[di];
      else for (int di = 0; di < (*sc).dim; di++) com[di] += (*top).fibers.e[ci].len * xi[di];
    }
  double wt = 1. / (*sc).top.v.len;
  for (int di = 0; di < (*sc).dim; di++) com[di] *= wt;
  for (int fi = 0; fi < (*sc).fxd.len; fi++)
    {
      int i = (*sc).fxd.e[fi];
      const double *xi = string_config_vertex_coords(sc, i);
      for (int di = 0; di < (*sc).dim; di++)
	{
	  double delx = xi[di] - com[di];
	  var += delx * delx;
	}
    }
  for (int mi = 0; mi < (*cmobile).len; mi++)
    {
      int ci = (*cmobile).e[mi];
      const double *xi = gsl_vector_const_ptr(c_data, (*c_map).e[mi]);
      for (int di = 0; di < (*sc).dim; di++)
	{
	  double delx = xi[di] - com[di];
	  delx *= delx;
	  if (top == NULL) {}
	  else delx *= (*top).fibers.e[ci].len;
      	  var += delx;
	}
    }
  return var / (*sc).top.v.len;
}

// NOTE: Make sure that 'c_map' and coordinate/force data (e.g. solvers) are updated concurrently (if necessary)
// anywhere this function is used. (Also, see below) This is more of a helper function for routines involving
// contracted/decontracted coordinate spaces than a stand-alone general purpose function (coordinate and force data
// aren't modified explicitly because it would be somewhat unwieldy to constantly allocate and deallocate vectors
// for each cluster.)
void sc_solver_expand_cluster_exp(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, int ci)
{
  // Define new clusters for each additional vertex in cluster 'ci'
  if ((*aux).fibers.e[ci].len > 1) {}
  else return;
  array_int old_fiber = (*aux).fibers.e[ci];
  // NOTE: this needs to be updated consistently with 'contr_nbrlist_expand' (e.g. whether 'ci' is removed or kept)
  for (int fi = 1; fi < old_fiber.len; fi++)
    {
      int mi = (*cmobile).len;
      add2array_int(cmobile, (*cmobile_map).len);
      add2array_int(cmobile_map, mi);
    }
  contr_nbrlist_expand(aux, ci, &((*sc).edge_wts));
}

void remove_embedding(array_int *emb, array_int *emb_inv, int emb_index)
{
  // List, List_addr, M: List_addr[List[mi]] = mi
  // Conjectured general approach: contract the embedding, then contract the base space
  //         - Remove mi from List
  //         - Set List_addr[ci] = -1, and (if mi < len(List)) List_addr[List[mi]] = mi
  //         - Remove ci from List_addr: this will impact List if len(M)-1 appears in List (and ci < len(M)-1)
  //         - If ci < len(M)-1, set List[List_addr[ci]] to ci.
  int inv_index = (*emb).e[emb_index];
  remove_array_int(emb, emb_index);
  if (emb_index < (*emb).len)
    {
      (*emb_inv).e[(*emb).e[emb_index]] = emb_index;
    }
  remove_array_int(emb_inv, inv_index);
  if (inv_index < (*emb_inv).len && (*emb_inv).e[inv_index] > -1) (*emb).e[(*emb_inv).e[inv_index]] = inv_index;  
}

void sc_minimizer_set_clusters_exp(string_config *sc, contr_nbrlist *aux, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double epsilon)
{
  double epsq = epsilon * epsilon;
  int mi = 0;
  int init_size = (*cmobile).len;
  while (mi < (*cmobile).len)
    {
      int ci = (*cmobile).e[mi];
      const double *xi = sc_solver_vertex_coords_exp2(sc, ci, &(*aux), &(*cmobile_map), &((*c_map)), x);
      int cii_ = -1;
      for (int ni = 0; ni < (*aux).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*aux).top.top.v.e[ci].e[ni];
	  if ((*cmobile_map).e[cii] < mi) {}
	  else continue;
	  const double *xii = sc_solver_vertex_coords_exp2(sc, cii, &(*aux), &(*cmobile_map), &((*c_map)), x);
	  double delsq = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delx = xi[di] - xii[di];
	      delsq += delx * delx;
	    }
	  if (delsq > epsq) {}
	  else
	    {
	      //printf("Found merging clusters: dist(%d,%d) = %g\n", ci, cii, sqrt(delsq));
	      cii_ = cii;
	      break;
	    }
	}
      if (cii_ > -1)
	{
	  //printf("Attempting to merge %d-%d\n", ci, cii_);
	  contr_nbrlist_merge_cluster(aux, ci, cii_);
	  remove_embedding(cmobile, cmobile_map, mi);
	  remove_array_int(c_map, mi);
	  //printf("(done)\n");
	}
      else mi += 1;
    }
}

// NOTE: this function is probably superfluous (because the 'rf' solver typically separates merging vertices)
void sc_minimizer_rf_reset(sc_minimizer_rf *scf, double epsilon, const gsl_multiroot_fsolver_type *gsl_fsolver_type)
{
  gsl_fsolver_type = gsl_fsolver_type == NULL ? gsl_multiroot_fsolver_hybrids : gsl_fsolver_type;
  string_config *sc = (*scf).sc;
  double epsq = epsilon * epsilon;
  sc_minimizer_set_clusters_exp((*scf).sc, &((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map), &((*scf).c_map), (*scf).solver->x, epsilon);
  // Re-read coordinates
  // Re-read coordinates and reset the solver
  (*scf).n_vars = (*scf).cmobile.len * (*sc).dim;
  gsl_vector_free((*scf).c_data);
  (*scf).c_data = gsl_vector_alloc((*scf).n_vars);
  read_mobile_coords_string_config(sc, &((*scf).top), &((*scf).cmobile), (*scf).c_data, &((*scf).c_map));
  gsl_multiroot_fsolver_free((*scf).solver);
  (*scf).solver = gsl_multiroot_fsolver_alloc(gsl_fsolver_type, (*scf).n_vars);
  (*scf).solver_data.n = (*scf).n_vars;
  gsl_multiroot_fsolver_set((*scf).solver, &((*scf).solver_data), (*scf).c_data);
}

int sc_minimizer_rf_iterate(sc_minimizer_rf *scf)
{
  sc_minimizer_rf_counter += 1;
  return gsl_multiroot_fsolver_iterate((*scf).solver);
}

int sc_minimizer_rf_relax(sc_minimizer_rf *scm, double tol)
{
  double tolsq = tol * tol;
  //(*scm).tlf_pars.core_radsq = tolsq / (*(*scm).sc).dim;
  int status_ = GSL_CONTINUE;
  int count = 0;
  while (count < MAX_ITER)
    {
      count += 1;
      int status = sc_minimizer_rf_iterate(scm);
      if (status == GSL_SUCCESS)
	{
	  gsl_vector *dx = gsl_multiroot_fsolver_dx((*scm).solver);
	  double *dx_ = gsl_vector_ptr(dx, 0);
	  double normsqdx = euclid_normsq(dx_, (*dx).size);
	  if (normsqdx > 0)
	    {
	      if (normsqdx < tolsq)
		{
		  status_ = GSL_SUCCESS;
		  break;
		}
	    }
	}
      else
	{
	  status_ = status;
	  break;
	}
    }
  return status_;
}

double aux_gsl_vector_dot(gsl_vector *a, gsl_vector *b)
{
  double dp = 0;
  double *a_ = gsl_vector_ptr(a, 0);
  double *b_ = gsl_vector_ptr(b, 0);
  int len = (*a).size;
  for (int i = 0; i < len; i++)
    {
      dp += a_[i] * b_[i];
    }
  return dp;
}

int sc_minimizer_rf_relax_diag(sc_minimizer_rf *scf, int n_steps, FILE *ofile)
{
  for (int i = 0; i < n_steps; i++)
    {
      int status = gsl_multiroot_fsolver_iterate((*scf).solver);
      fprintf(ofile, "%d ", status);
      for (int ci = 0; ci < (*(*scf).solver).x->size; ci++) fprintf(ofile, "%g ", gsl_vector_get((*(*scf).solver).x, ci));
      double total_length = total_length_func((*(*scf).solver).x, (*scf).pars);
      gsl_vector *f = gsl_multiroot_fsolver_f((*scf).solver);
      double norm_f = sqrt(aux_gsl_vector_dot(f, f));
      fprintf(ofile, "%g %g %d %d %d\n", total_length, norm_f, total_length_func_counter, grad_total_length_func_counter, comp_total_length_func_counter);
    }
}

void sc_minimizer_mm_init_exp(sc_minimizer_mm *scm, string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *x, double stepsize, double tol)
{
  (*scm).sc = sc;
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  prep_contr_list(&((*sc).mobile), &((*scm).cmobile), &((*scm).cmobile_map), (*sc).top.v.len);
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  total_length_func_pars *tlf_pars = (total_length_func_pars *) malloc(sizeof(total_length_func_pars));
  if ((*scm).cmobile.len > 0)
  {
	  (*scm).c_map.len = (*scm).cmobile.len;
	  int n_vars = (*scm).cmobile.len * (*sc).dim;
	  (*scm).c_data = gsl_vector_alloc(n_vars);
	  sc_solver_read_mobile_coords_exp2(sc, top, cmobile, cmobile_map, c_map, x, &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), &((*scm).c_map), (*scm).c_data);
	  total_length_func_pars_init(tlf_pars, sc, &((*scm).c_map), stepsize);
	  (*scm).pars = tlf_pars;
	  (*tlf_pars).top = &((*scm).top);
	  (*tlf_pars).cmobile = &((*scm).cmobile);
	  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
	  (*tlf_pars).core_radsq = 1e-32;
	  (*scm).n_vars = n_vars;
	  //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
	  //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, n_vars);
	  (*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
	  (*scm).solver_data.f = total_length_func;
	  (*scm).solver_data.df = grad_total_length_func;
	  (*scm).solver_data.fdf = comp_total_length_func;
	  (*scm).solver_data.params = (*scm).pars;
	  (*scm).solver_data.n = n_vars;
	  gsl_multimin_fdfminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, stepsize, tol);
  }
  else
  {
	  int n_vars = 0;
	  (*scm).c_data = NULL;
	  // This may be unnecessary
	  total_length_func_pars_init(tlf_pars, sc, &((*scm).c_map), stepsize);
	  (*tlf_pars).top = &((*scm).top);
	  (*tlf_pars).cmobile = &((*scm).cmobile);
	  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
	  (*tlf_pars).core_radsq = 1e-32;
	  (*scm).pars = tlf_pars;
	  (*scm).n_vars = 0;
	  //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
	  //(*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, n_vars);
	  (*scm).solver = NULL;
  }
}

void sc_minimizer_mm_init(sc_minimizer_mm *scm, string_config *sc)
{
  double stepsize = sqrt(string_config_var_x(sc)) * 0.01;
  double tol = 1e-2;
  (*scm).sc = sc;
  contr_nbrlist_init(&((*scm).top), &((*sc).top), &((*sc).edge_wts));
  prep_contr_list(&((*sc).mobile), &((*scm).cmobile), &((*scm).cmobile_map), (*sc).top.v.len);
  contract_leaves(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  contract_elbows(&((*scm).top), NULL, &((*scm).cmobile), &((*scm).cmobile_map));
  array_int_init(&((*scm).c_map), (*scm).cmobile.len);
  (*scm).c_map.len = (*scm).cmobile.len;
  total_length_func_pars *tlf_pars = (total_length_func_pars *) calloc(1, sizeof(total_length_func_pars));
  (*scm).pars = tlf_pars;
  if ((*scm).cmobile.len > 0)
  {
	  int n_vars = (*scm).cmobile.len * (*sc).dim;
	  (*scm).n_vars = n_vars;
	  (*scm).c_data = gsl_vector_alloc(n_vars);
	  (*scm).solver_data.f = total_length_func;
	  (*scm).solver_data.df = grad_total_length_func;
	  (*scm).solver_data.fdf = comp_total_length_func;
	  (*scm).solver_data.params = (*scm).pars;
	  (*scm).solver_data.n = n_vars;
	  read_mobile_coords_string_config(sc, &((*scm).top), &((*scm).cmobile), (*scm).c_data, &((*scm).c_map));
	  double *c_data_ = gsl_vector_ptr((*scm).c_data, 0);
	  for (int di = 0; di < (*(*scm).c_data).size; di++)
	    {
	      c_data_[di] += stepsize * (rnd() - 0.5);
	    }
	  total_length_func_pars_init(tlf_pars, sc, &((*scm).c_map), stepsize);
	  (*tlf_pars).top = &((*scm).top);
	  (*tlf_pars).cmobile = &((*scm).cmobile);
	  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
	  (*tlf_pars).core_radsq = 1e-32;
	  (*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
	  gsl_multimin_fdfminimizer_set((*scm).solver, &((*scm).solver_data), (*scm).c_data, stepsize, tol);
  }
  else
  {
	  int n_vars = 0;
	  (*scm).n_vars = n_vars;
	  (*scm).c_data = NULL;
	  total_length_func_pars_init(tlf_pars, sc, &((*scm).c_map), stepsize);
	  (*tlf_pars).top = &((*scm).top);
	  (*tlf_pars).cmobile = &((*scm).cmobile);
	  (*tlf_pars).cmobile_map = &((*scm).cmobile_map);
	  (*tlf_pars).core_radsq = 1e-32;
	  (*scm).solver = NULL;
  }
}

void free_sc_minimizer_mm(sc_minimizer_mm *scm)
{
  free((*scm).pars);
  free_array_int(&((*scm).c_map));
  if ((*scm).solver != NULL) gsl_multimin_fdfminimizer_free((*scm).solver);
  if ((*scm).c_data != NULL) gsl_vector_free((*scm).c_data);
  free_array_int(&((*scm).cmobile));
  free_array_int(&((*scm).cmobile_map));
  free_contr_nbrlist(&((*scm).top));
}

void sc_minimizer_mm_contract(sc_minimizer_mm *scf, double epsilon, const gsl_multimin_fdfminimizer_type *gsl_alg_type, double stepsize, double tol)
{
  gsl_alg_type = gsl_alg_type == NULL ? gsl_minimizer_type_mm : gsl_alg_type;
  string_config *sc = (*scf).sc;
  double epsq = epsilon * epsilon;
  sc_minimizer_set_clusters_exp((*scf).sc, &((*scf).top), &((*scf).cmobile), &((*scf).cmobile_map), &((*scf).c_map), (*scf).solver->x, epsilon);
  if ((*scf).cmobile.len > 0)
  {
	  (*scf).n_vars = (*scf).cmobile.len * (*sc).dim;
	  gsl_vector *nc_data = gsl_vector_alloc((*scf).n_vars);
	  int base = 0;
	  for (int cmi = 0; cmi < (*scf).cmobile.len; cmi++)
	    {
	      int ci = (*scf).cmobile.e[cmi];
	      double *xi0 = gsl_vector_ptr((*scf).solver->x, (*scf).c_map.e[cmi]);
	      double *xi = gsl_vector_ptr(nc_data, base);
	      for (int di = 0; di < (*sc).dim; di++)
		{
		  xi[di] = xi0[di];
		}
	      (*scf).c_map.e[cmi] = base;
	      base += (*sc).dim;
	    }
	  gsl_vector_free((*scf).c_data);
	  (*scf).c_data = nc_data;
	  gsl_multimin_fdfminimizer_free((*scf).solver);
	  (*scf).solver = gsl_multimin_fdfminimizer_alloc(gsl_alg_type, (*scf).n_vars);
	  (*scf).solver_data.n = (*scf).n_vars;
	  gsl_multimin_fdfminimizer_set((*scf).solver, &((*scf).solver_data), (*scf).c_data, stepsize, tol);
  }
  else
  {
	  (*scf).n_vars = 0;
	  gsl_vector_free((*scf).c_data);
	  (*scf).c_data = NULL;
	  gsl_multimin_fdfminimizer_free((*scf).solver);
	  (*scf).solver = NULL;
  }
}

int sc_minimizer_test_stability_mm(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol)
{
  sc_minimizer_mm scm;
  sc_minimizer_mm_init_exp(&scm, sc, top, cmobile, cmobile_addr, c_map, x, 1e-4, 1e-2);
  gsl_vector *x0 = gsl_vector_alloc((*x).size);
  gsl_vector_memcpy(x0, x);
  sc_minimizer_mm_relax(&scm, 0.1 * tol);
  // Measure the displacement from x0
  // double sc_solver_distsq_exp(string_config *sc, contr_nbrlist *top0, array_int *cmobile0, array_int *cmobile_map0, array_int *c_map0, const gsl_vector *x0, contr_nbrlist *top1, array_int *cmobile1, array_int *cmobile_map1, array_int *c_map1, const gsl_Vector *x1)
  double dispsq = sc_solver_distsq_exp(sc, top, cmobile, cmobile_addr, c_map, x, &(scm.top), &(scm.cmobile), &(scm.cmobile_map), &(scm.c_map), scm.solver->x);
  gsl_vector_free(x0);
  free_sc_minimizer_mm(&scm);
  return dispsq > (tol * tol);
}

// Alternatively: consider measuring the expansion of individual clusters
int sc_minimizer_test_stability_sd(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_addr, array_int *c_map, const gsl_vector *x, double tol)
{
  int n_vars = (*sc).dim * (*sc).mobile.len;
  int n_steps = 100;
  int max_total_wt = string_config_max_total_weight(sc);
  n_steps *= max_total_wt;
  gsl_vector *exp_data = gsl_vector_alloc(n_vars);
  // Read initial coordinates (in uncontracted space)
  int base = 0;
  array_int aux_cmap;
  array_int_init(&aux_cmap, (*sc).mobile.len);
  aux_cmap.len = (*sc).mobile.len;
  for (int mi = 0; mi < (*sc).mobile.len; mi++)
    {
      int i = (*sc).mobile.e[mi];
      int ci = (*top).map.e[i];
      int cmi = (*cmobile_addr).e[ci];
      const double *xci;
      double *xi = gsl_vector_ptr(exp_data, base);
      if (cmi > -1) xci = gsl_vector_const_ptr(x, (*c_map).e[cmi]);
      else xci = string_config_vertex_coords(sc, i);
      for (int di = 0; di < (*sc).dim; di++) xi[di] = xci[di];
      aux_cmap.e[mi] = base;
      base += (*sc).dim;
    }
  // Attempt to relax configuration over n_steps steps
  gsl_vector *x0 = gsl_vector_alloc(n_vars);
  gsl_vector_memcpy(x0, exp_data);
  gsl_vector *f = gsl_vector_alloc(n_vars);
  for (int si = 0; si < n_steps; si++)
    {
      string_config_compute_force_mobile_coords(sc, &aux_cmap, exp_data, f);
      gsl_blas_daxpy(tol, f, exp_data);
    }
  // Check net displacement
  gsl_vector_sub(exp_data, x0);
  double *dx = gsl_vector_ptr(exp_data, 0);
  double norm_dispsq = euclid_normsq(dx, (*exp_data).size); // Consider normalizing this by the number of cluster vertices
  double thr = tol * (n_steps * tol + max_total_wt); // NOTE: consider adjusting this: replace the constant with an empirically determined envelope (as a function of the 'time step' and local graph properties near singular configurations)
  thr *= thr;
  int status = norm_dispsq > thr;
  // Free auxiliary data
  gsl_vector_free(exp_data);
  gsl_vector_free(x0);
  gsl_vector_free(f);
  free_array_int(&aux_cmap);
  return status;
}

void sc_minimizer_mm_reset(sc_minimizer_mm *scm)
{
  string_config *sc = (*scm).sc;
  int n_vars = (*sc).dim * (*sc).mobile.len;
  if (n_vars != (*scm).n_vars)
    {
      gsl_vector *nc_data = gsl_vector_alloc(n_vars);
      sc_solver_expand(sc, &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), &((*scm).c_map), (*scm).solver->x, nc_data); 
      gsl_multimin_fdfminimizer_free((*scm).solver);
      (*scm).n_vars = n_vars;
      (*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, (*scm).n_vars);
      (*scm).solver_data.n = n_vars;
      double init_step = sqrt(sc_minimizer_var_mobile_coords(sc, &((*scm).top), &(scm->cmobile), &(scm->c_map), scm->solver->x)) * 0.01;
      gsl_multimin_fdfminimizer_set((*scm).solver, &((*scm).solver_data), nc_data, init_step, 0.01);
      gsl_vector_free((*scm).c_data);
      (*scm).c_data = nc_data;
    }
}

int sc_minimizer_mm_solve(sc_minimizer_mm *scm, double tol)
{
  gsl_vector *x0 = gsl_vector_alloc((*scm).n_vars);
  gsl_vector_memcpy(x0, (*scm).solver->x);
  contr_nbrlist t0;
  transcribe_contr_nbrlist(&((*scm).top), &t0); 
  array_int cm0;
  transcribe_array_int(&((*scm).cmobile), &cm0);
  array_int cmm0;
  transcribe_array_int(&((*scm).cmobile_map), &cmm0);
  array_int cmap0;
  transcribe_array_int(&((*scm).c_map), &cmap0);
  double tolsq = tol * tol;
  int rstatus = GSL_CONTINUE;
  int count = 0;
  while (count < MAX_ITER)
    {
      count += 1;
      int status = sc_minimizer_mm_relax(scm, tol);
      // RESUME: check that status is acceptable before progressing (i.e. that the 'relaxation' hasn't terminated in a failure state)
     
      char realloc_flag = 0;
      if (contr_nbrlist_fiber_corresp(&t0, &((*scm).top), &((*(*scm).sc).mobile))) {}
      else realloc_flag = 1;
      double dispsq = sc_solver_distsq_exp((*scm).sc, &t0, &cm0, &cmm0, &cmap0, x0, &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), &((*scm).c_map), (*scm).solver->x);
      if (dispsq > tolsq && (status != GSL_ENOPROG || realloc_flag == 1)) {}
      else
	{
	  rstatus = status;
	  break;
	}
      if (realloc_flag)
	{
	  // NOTE: there's probably a slightly more efficient approach to this (e.g. 'aligning' contraction t0 with (*scm).top)
	  // but this isn't (or shouldn't be) the rate-limiting step.
	  free_contr_nbrlist(&t0);
	  free_array_int(&cm0);
	  free_array_int(&cmm0);
	  free_array_int(&cmap0);
	  gsl_vector_free(x0);
	  transcribe_contr_nbrlist(&((*scm).top), &t0);
	  transcribe_array_int(&((*scm).cmobile), &cm0);
	  transcribe_array_int(&((*scm).cmobile_map), &cmm0);
	  transcribe_array_int(&((*scm).c_map), &cmap0);
	  x0 = gsl_vector_alloc((*scm).n_vars);
	}
      gsl_vector_memcpy(x0, (*scm).c_data);
      // RESET SOLVER AT CURRENT POSITION
      sc_minimizer_mm_reset(scm);
      // Consider perturbing mobile coordinates by random vector.
      for (int i = 0; i < 3; i++) sc_minimizer_mm_iterate(scm);
    }
  // Free auxiliary structures
   free_contr_nbrlist(&t0);
   free_array_int(&cm0);
   free_array_int(&cmm0);
   free_array_int(&cmap0);
   gsl_vector_free(x0);  
   return rstatus;
}

int sc_minimizer_mm_relax(sc_minimizer_mm *scm, double tol)
{
  if ((*scm).n_vars > 0) {}
  else
    {
      printf("done (trivial case of sc_minimizer_mm)\n");
      return GSL_SUCCESS;
    }
  string_config *sc = (*scm).sc;
  double tolsq = tol * tol;
  double f0 = gsl_multimin_fdfminimizer_minimum((*scm).solver);
  int status_ = GSL_CONTINUE;
  int count = 0;
  double merge_rad = 10 * tol;
  while (count < MAX_ITER)
    {
      count += 1;
      int status;
      for (int i = 0; i < 5; i++) status = sc_minimizer_mm_iterate(scm);
      //double f1 = gsl_multimin_fdfminimizer_minimum((*scm).solver);
      //gsl_vector *dx = gsl_multimin_fdfminimizer_dx((*scm).solver);
      //double *dx_ = gsl_vector_ptr(dx, 0);
      double *df_ = gsl_vector_ptr(gsl_multimin_fdfminimizer_gradient((*scm).solver), 0);
      //double df = f1 - f0;
      double normsqdx = euclid_normsq(df_, (*scm).n_vars);
      //if (normsqdx > 0)
      //{
      if (normsqdx < tolsq)
	{
	  status_ = GSL_SUCCESS;
	  break;
	}
      //}
      //f0 = f1;
      // Check for merging clusters (and re-allocate the solver if necessary)
      int init_len = (*scm).n_vars;
      sc_minimizer_mm_check_merging(scm, merge_rad);
      if ((*scm).solver == NULL)
	{
	  printf("Solver set to NULL (minimization terminated)\n");
	  status_ = GSL_ENOPROG;
	  break;
	}
      if (status == GSL_SUCCESS && (*scm).solver != NULL) {}
      else
	{
	  //printf("(mm_relax): Status = %d after %d iterations\n", status, count);
	  status_ = status;
	  if (init_len == (*scm).n_vars) break;
	}
    }
  return status_;
}

void sc_minimizer_check_merging(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector **c_data, double merge_rad)
{
  int nc0 = (*top).top.top.v.len;
  sc_minimizer_set_clusters_exp(sc, top, cmobile, cmobile_map, c_map, *c_data, merge_rad);
  int nc1 = (*top).top.top.v.len;
  if (nc0 != nc1)
    {
      // Reallocate/condense coordinate data and update c_map
      int n_vars = (*cmobile).len * (*sc).dim;
      gsl_vector *nc_data = gsl_vector_alloc(n_vars);
      int base = 0;
      for (int cmi = 0; cmi < (*cmobile).len; cmi++)
	{
	  double *xi = gsl_vector_ptr(nc_data, base);
	  double *xi0 = gsl_vector_ptr(*c_data, (*c_map).e[cmi]);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      xi[di] = xi0[di];
	    }
	  (*c_map).e[cmi] = base;
	  base += (*sc).dim;
	}
      gsl_vector_free(*c_data);
      (*c_data) = nc_data;
    }
}

void sc_minimizer_mm_check_merging(sc_minimizer_mm *scm, double merge_rad)
{
  if ((*scm).n_vars > 0) {}
  else return;
  string_config *sc = (*scm).sc;
  contr_nbrlist *top = &((*scm).top);
  array_int *cmobile_map = &((*scm).cmobile_map);
  array_int *cmobile = &((*scm).cmobile);
  array_int *c_map = &((*scm).c_map);
  gsl_vector *c_data = gsl_vector_alloc((*scm).n_vars);
  gsl_vector_memcpy(c_data, (*scm).solver->x);
  int init_n_vs = (*cmobile).len;
  sc_minimizer_check_merging(sc, top, cmobile, cmobile_map, c_map, &c_data, merge_rad);
  int final_n_vs = (*cmobile).len;
  if (final_n_vs != init_n_vs)
    {
      gsl_multimin_fdfminimizer_free((*scm).solver);
      gsl_vector_free((*scm).c_data);
      if (final_n_vs > 0)
	{
	  // Reset the solver
	  (*scm).n_vars = (*c_data).size;
	  (*scm).solver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, (*scm).n_vars);
	  double stepsize = sqrt(sc_minimizer_var_mobile_coords((*scm).sc, top, cmobile, c_map, c_data)) / (*cmobile).len * 0.01;
	  double tol = 1e-2;
	  (*scm).solver_data.n = (*c_data).size;
	  (*scm).c_data = c_data;
	  gsl_multimin_fdfminimizer_set((*scm).solver, &((*scm).solver_data), c_data, stepsize, tol);
	  double *df = gsl_vector_ptr(gsl_multimin_fdfminimizer_gradient((*scm).solver), 0);
	  double norm_gradient = sqrt(euclid_normsq(df, (*scm).n_vars));
	}
      else
	{
	  (*scm).n_vars = 0;
	  (*scm).solver = NULL;
	  (*scm).c_data = NULL;
	  gsl_vector_free(c_data);
	}
    }
  else gsl_vector_free(c_data);
}

int sc_minimizer_mm_iterate(sc_minimizer_mm *scm)
{
  if ((*scm).n_vars > 0)
    {
      sc_minimizer_mm_counter += 1;
      return gsl_multimin_fdfminimizer_iterate((*scm).solver);
    }
  else return GSL_ENOPROG;
}

int sc_minimizer_mm_relax_diag(sc_minimizer_mm *scm, int n_steps, FILE *ofile)
{
  if ((*scm).n_vars > 0) {}
  else
  {
	printf("Error: sc_minimizer_mm instance has no mobile coordinates\n");
	return GSL_SUCCESS;
  }
  total_length_func_pars tlf_pars;
  tlf_pars.top = &((*scm).top);
  tlf_pars.sc = (*scm).sc;
  tlf_pars.cmobile = &((*scm).cmobile);
  tlf_pars.cmobile_map = &((*scm).cmobile_map);
  tlf_pars.c_map = &((*scm).c_map);
  tlf_pars.core_radsq = 1e-32;
  for (int i = 0; i < n_steps; i++)
    {
      int status = sc_minimizer_mm_iterate(scm);
      gsl_vector *x = (*(*scm).solver).x;
      fprintf(ofile, "%d ", status);
      for (int ci = 0; ci < (*x).size; ci++)
	{
	  fprintf(ofile, "%g ", gsl_vector_get(x, ci));
	}
      gsl_vector *g = gsl_multimin_fdfminimizer_gradient(scm->solver);
      double norm_g = sqrt(aux_gsl_vector_dot(g, g));
      fprintf(ofile, "%g %g %d %d %d\n", norm_g, total_length_func(x, &tlf_pars), total_length_func_counter, grad_total_length_func_counter, comp_total_length_func_counter);
    }
  return GSL_SUCCESS;
}
 
// <<< sc_Lagrange_solver >>

int Lagrange_solver_f(const gsl_vector *c_data, void *params, gsl_vector *f)
{
  // The gradient of the conventional energy function:
  // 	E = tau * (length - sc.L0) - ext_f . x_ef
  // dE = dtau * (length - sc.L0) + tau * dlength - ext_f . dx_ef
  //
  sc_Lagrange_solver *scs = (sc_Lagrange_solver *) params;
  string_config *sc = (*scs).sc;
  
  double len_disc = Lagrange_solver_f_exp((*scs).sc, &((*scs).top), &((*scs).cmobile), &((*scs).cmobile_map), &((*scs).c_map), c_data, f);
  //comp_total_length_func_exp(c_data, sc, &((*scs).c_map), &len_disc, f);
  gsl_vector_scale(f, gsl_vector_get(c_data, (*scs).tau_addr));
  len_disc -= (*scs).len_disc;
  gsl_vector_set(f, (*scs).tau_addr, len_disc);
  int ci = (*scs).top.map.e[(*scs).ext_f_site];
  int mci = (*scs).cmobile_map.e[ci];
  if (mci == -1)
    {
      printf("Error: this seems like a logical impossibility... or at least that it *should* be a logical impossibility. (Bye.)\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      double *fi = gsl_vector_ptr(f, (*scs).c_map.e[mci]);
      for (int di = 0; di < (*sc).dim; di++) fi[di] -= (*scs).ext_f[di];
    }
  return GSL_SUCCESS;
}

void sub_lower_triangular(gsl_matrix *dest, gsl_matrix *src)
{
  for (int i = 0; i < (*dest).size1; i++)
    {
      for (int j = 0; j <= i; j++)
	{
	  double a = gsl_matrix_get(dest, i, j) - gsl_matrix_get(src, i, j);
	  gsl_matrix_set(dest, i, j, a);
	}
    }
}

void add_lower_triangular(gsl_matrix *dest, gsl_matrix *src)
{
  for (int i = 0; i < (*dest).size1; i++)
    {
      for (int j = 0; j <= i; j++)
	{
	  double a = gsl_matrix_get(dest, i, j) + gsl_matrix_get(src, i, j);
	  gsl_matrix_set(dest, i, j, a);
	}
    }
}

void Lagrange_solver_df_compute_incr(gsl_matrix *incr, const double *xi, const double *xii, double *delx, int wt, double *lensq, double *delta_f_i, double tau)
{
  int dim = (*incr).size1;
  //double delx[dim];
  //double delxsq = 0;
  double dist;
  double inv_delxsq = 1. / (*lensq); 
  double wt_inv_dist = wt * sqrt(inv_delxsq);
  double diag_coeff = wt_inv_dist * tau;
  gsl_matrix_set_zero(incr);
  if (dim < 5)
    {
      double op_coeff = inv_delxsq * diag_coeff;
      for (int di = 0; di < dim; di++)
	{
	  delta_f_i[di] = delx[di] * wt_inv_dist;
	  gsl_matrix_set(incr, di, di, diag_coeff - op_coeff * delx[di] * delx[di]);
	  for (int dii = 0; dii < di; dii++) 
	    {
	      double a = -op_coeff * delx[di] * delx[dii];
	      gsl_matrix_set(incr, di, dii, a);
	    }
	}
    }
  else
    {
      double C1 = sqrt(diag_coeff * inv_delxsq);
      for (int di = 0; di < dim; di++) 
	{
	  delta_f_i[di] = wt_inv_dist * delx[di];
	  delx[di] *= C1;
	}
      for (int di = 0; di < dim; di++)
	{
	  gsl_matrix_set(incr, di, di, diag_coeff - delx[di] * delx[di]);
	  for (int dii = 0; dii < di; dii++) 
	    {
	      double a = -delx[di] * delx[dii];
	      gsl_matrix_set(incr, di, dii, a);
	    }
	}
    }
}

int Lagrange_solver_df(const gsl_vector *c_data, void *params, gsl_matrix *df)
{
  sc_Lagrange_solver *scs = (sc_Lagrange_solver *) params;
  string_config *sc = (*scs).sc;
  gsl_matrix_set_zero(df);
  gsl_vector_view dxdtau = gsl_matrix_row(df, (*scs).tau_addr);
  (*scs).mobile_length = 0;
  for (int mi = 0; mi < (*scs).cmobile.len; mi++)
    {
      int i = (*scs).cmobile.e[mi];
      double *f_i = gsl_vector_ptr(&(dxdtau.vector), (*scs).c_map.e[mi]);
      const double *x_i = gsl_vector_const_ptr(c_data, (*scs).c_map.e[mi]);
      gsl_matrix_view block_i = gsl_matrix_submatrix(df, (*scs).c_map.e[mi], 0, (*sc).dim, (*scs).tau_addr);
      gsl_matrix_view block_i_i = gsl_matrix_submatrix(&(block_i.matrix), 0, (*scs).c_map.e[mi], (*sc).dim, (*sc).dim);
      for (int ni = 0; ni < (*scs).top.top.top.v.e[i].len; ni++)
	{
	  int ii = (*scs).top.top.top.v.e[i].e[ni];
	  int mii = (*scs).cmobile_map.e[ii];
	  if (mii < mi)
	    {
	      const double *x_ii;
	      if (mii == -1) x_ii = string_config_vertex_coords(sc, (*scs).top.fibers.e[ii].e[0]);
	      else x_ii = gsl_vector_const_ptr(c_data, (*scs).c_map.e[mii]);
	      double disp[(*sc).dim];
	      double distsq = 0;
	      for (int di = 0; di < (*sc).dim; di++)
		{
		  disp[di] = x_ii[di] - x_i[di];
		  distsq += disp[di] * disp[di];
		}
	      if (distsq > 0) {}
	      else continue;
	      gsl_matrix *incr = gsl_matrix_alloc((*sc).dim, (*sc).dim);
	      double edge_len_sq = distsq;
	      double delta_f_i[(*sc).dim];
	      Lagrange_solver_df_compute_incr(incr, x_i, x_ii, &disp[0], (*scs).top.top.edge_wts.e[i].e[ni], &edge_len_sq, &delta_f_i[0], gsl_vector_get(c_data, (*scs).tau_addr)); 
	      (*scs).mobile_length += (*scs).top.top.edge_wts.e[i].e[ni] * sqrt(edge_len_sq);
	      for (int di = 0; di < (*sc).dim; di++) f_i[di] -= delta_f_i[di];
	      add_lower_triangular(&(block_i_i.matrix), incr);
	      if (mii == -1) {}
	      else
		{
		  double *f_ii = gsl_vector_ptr(&(dxdtau.vector), (*scs).c_map.e[mii]);
		  for (int di = 0; di < (*sc).dim; di++) f_ii[di] += delta_f_i[di];
		  gsl_matrix_view block_ii_ii = gsl_matrix_submatrix(df, (*scs).c_map.e[mii], (*scs).c_map.e[mii], (*sc).dim, (*sc).dim);
		  gsl_matrix_view block_i_ii = gsl_matrix_submatrix(&(block_i.matrix), 0, (*scs).c_map.e[mii], (*sc).dim, (*sc).dim);
		  add_lower_triangular(&(block_ii_ii.matrix), incr);
		  sub_lower_triangular(&(block_i_ii.matrix), incr);
		}
	      gsl_matrix_free(incr);
	    }
	}
    }
  // Set symmetric terms in df
  for (int mi = 0; mi < (*scs).cmobile.len; mi++)
    {
      for (int mii = 0; mii < mi; mii++)
	{
	  gsl_matrix_view block_i_ii = gsl_matrix_submatrix(df, (*scs).c_map.e[mi], (*scs).c_map.e[mii], (*sc).dim, (*sc).dim);
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      for (int dii = 0; dii < di; dii++)
		{
		  double a = gsl_matrix_get(&(block_i_ii.matrix), di, dii);
		  gsl_matrix_set(&(block_i_ii.matrix), dii, di, a);
		}
	    }
	}
    }
  for (int i = 0; i < (*df).size1; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  double a = gsl_matrix_get(df, i, j);
	  gsl_matrix_set(df, j, i, a);
	}
    }
  return GSL_SUCCESS;
}

int Lagrange_solver_fdf(const gsl_vector *c_data, void *params, gsl_vector *f, gsl_matrix *df)
{
  int status = Lagrange_solver_df(c_data, params, df);
  if (status == GSL_SUCCESS) {}
  else return status;
  // Transcribe the last row of df into f and scale by tau
  sc_Lagrange_solver *scs = (sc_Lagrange_solver *) params;
  gsl_vector_view lr = gsl_matrix_subrow(df, (*scs).tau_addr, 0, (*scs).n_vars);
  gsl_vector_memcpy(f, &(lr.vector));
  gsl_vector_scale(f, gsl_vector_get(c_data, (*scs).tau_addr));
  int ci = (*scs).top.map.e[(*scs).ext_f_site];
  int mci = (*scs).cmobile_map.e[ci];
  if (mci > -1)
    {
      double *f_i = gsl_vector_ptr(f, (*scs).c_map.e[mci]);
      for (int di = 0; di < (*(*scs).sc).dim; di++) f_i[di] -= (*scs).ext_f[di];
    }
  else
    {
      printf("Something weird happened! a;s;oja;s;a;s\n");
      exit(EXIT_FAILURE);
    }
  gsl_vector_set(f, (*scs).tau_addr, (*scs).mobile_length - (*scs).len_disc);
  return GSL_SUCCESS;
}

// The following function is mainly intended for debugging
double sc_solver_compute_total_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data)
{
  double len = 0;
  for (int ci = 0; ci < (*top).top.top.v.len; ci++)
    {
      const double *xci;
      int mci = (*cmobile_map).e[ci];
      if (mci > -1) xci = gsl_vector_const_ptr(c_data, (*c_map).e[mci]);
      else
	{
	  int ri = (*top).fibers.e[ci].e[0];
	  xci = string_config_vertex_coords(sc, ri);
	}
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  if (cii < ci) continue;
	  const double *xcii;
	  int mcii = (*cmobile_map).e[cii];
	  if (mcii > -1) xcii = gsl_vector_const_ptr(c_data, (*c_map).e[mcii]);
	  else xcii = string_config_vertex_coords(sc, (*top).fibers.e[cii].e[0]);
	  double incr = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delx = xcii[di] - xci[di];
	      incr += delx * delx;
	    }
	  len += (*top).top.edge_wts.e[ci].e[ni] * sqrt(incr);
	}
    }
  return len;
}

double sc_solver_compute_mobile_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data)
{
  double len = 0;
  for (int mi = 0; mi < (*cmobile).len; mi++)
    {
      int ci = (*cmobile).e[mi];
      const double *xi = gsl_vector_const_ptr(c_data, (*c_map).e[mi]);
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  int mcii = (*cmobile_map).e[cii];
	  if (mcii < mi) {}
	  else continue;
	  const double *xii;
	  if (mcii > -1)
	    {
	      xii = gsl_vector_const_ptr(c_data, (*c_map).e[mcii]);
	    }
	  else xii = string_config_vertex_coords(sc, (*top).fibers.e[cii].e[0]);
	  double incr = 0;
	  double delx[(*sc).dim];
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] = xi[di] - xii[di];
	      incr += delx[di] * delx[di];
	    }
	  incr = sqrt(incr);
	  len += (*top).top.edge_wts.e[ci].e[ni] * incr;
	}
    }
  return len;
}

double sc_solver_compute_fixed_length_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data)
{
  double len = 0;
  for (int mi = 0; mi < (*sc).fxd.len; mi++)
    {
      int i = (*sc).fxd.e[mi];
      const double *xi = string_config_vertex_coords(sc, i);
      int ci = (*top).map.e[i];
      if (i == (*top).fibers.e[ci].e[0]) {}
      else
	{
	  int ri = (*top).fibers.e[ci].e[0];
	  printf("Inconsistency found in contr_nbrlist! c.f. %d, %d\n", i, ri);
	  exit(EXIT_FAILURE);
	}
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  int mcii = (*cmobile_map).e[cii];
	  if ((mcii == -1) && (cii < ci)) {}
	  else continue;
	  int rii = (*top).fibers.e[cii].e[0];
	  const double *xii = string_config_vertex_coords(sc, rii);
	  double incr = 0;
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      double delx = xi[di] - xii[di];
	      incr += delx * delx;
	    }
	  incr = (*top).top.edge_wts.e[ci].e[ni] * sqrt(incr);
	  len += incr;
	}
    }
  return len;
}

double Lagrange_solver_f_exp(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, const gsl_vector *c_data, gsl_vector *f)
{
  gsl_vector_set_zero(f);
  double len = 0;
  for (int mi = 0; mi < (*cmobile).len; mi++)
    {
      int ci = (*cmobile).e[mi];
      const double *xi = gsl_vector_const_ptr(c_data, (*c_map).e[mi]);
      double *fi = gsl_vector_ptr(f, (*c_map).e[mi]);
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  int mcii = (*cmobile_map).e[cii];
	  if (mcii < mi) {}
	  else continue;
	  const double *xii;
	  double *fii = NULL;
	  if (mcii > -1)
	    {
	      xii = gsl_vector_const_ptr(c_data, (*c_map).e[mcii]);
	      fii = gsl_vector_ptr(f, (*c_map).e[mcii]);
	    }
	  else xii = string_config_vertex_coords(sc, (*top).fibers.e[cii].e[0]);
	  double incr = 0;
	  double delx[(*sc).dim];
	  for (int di = 0; di < (*sc).dim; di++)
	    {
	      delx[di] = xi[di] - xii[di];
	      incr += delx[di] * delx[di];
	    }
	  incr = sqrt(incr);
	  len += (*top).top.edge_wts.e[ci].e[ni] * incr;
	  if (incr > 0)
	    {
	      incr = (*top).top.edge_wts.e[ci].e[ni] / incr;
	      for (int di = 0; di < (*sc).dim; di++) delx[di] *= incr;
	      for (int di = 0; di < (*sc).dim; di++) fi[di] += delx[di];
	      if (fii != NULL)
		{
		  for (int di = 0; di < (*sc).dim; di++) fii[di] -= delx[di];
		}
	    }
	}
    }
  return len;
}

void sc_Lagrange_solver_init(sc_Lagrange_solver *scs, string_config *sc, int ext_f_site, double *ext_f, double L)
{
  sc_minimizer_mm scm;
  sc_minimizer_mm_init_heuristic(&scm, sc, ext_f_site, ext_f, L, 1.0);
  heuristic_pars *hpars = (heuristic_pars *) scm.pars;
  sc_minimizer_mm_solve(&scm, 1e-8); // RESUME: determine an acceptable accuracy/tolerance/precision
  double stepsize = 0.002 * sqrt(string_config_var_x(sc));
  for (int i = 0; i < 2; i++)
    {
      (*hpars).k *= 4;
      gsl_vector_memcpy(scm.c_data, scm.solver->x);
      gsl_multimin_fdfminimizer_set(scm.solver, &(scm.solver_data), scm.c_data, stepsize, 1e-2);
      sc_minimizer_mm_solve(&scm, 1e-8);
      stepsize *= 0.2;
    }
  sc_Lagrange_solver_from_minimizer(scs, &scm, ext_f_site, ext_f, L);
  free_sc_minimizer_mm(&scm);
}

void sc_Lagrange_solver_from_minimizer(sc_Lagrange_solver *scs, sc_minimizer_mm *scm, int ext_f_site, double *ext_f, double L)
{
  string_config *sc = (*scm).sc;
  (*scs).sc = sc;
  (*scs).ext_f_site = ext_f_site;
  (*scs).ext_f = ext_f;
  (*scs).L = L;
  // Consider initializing the solver with minimal total length
  (*scs).solver_type = POWELL_SOLVER;
  (*scs).tau_addr = (*scm).n_vars;
  (*scs).n_vars = (*scs).tau_addr + 1;
  transcribe_contr_nbrlist(&((*scm).top), &((*scs).top));
  transcribe_array_int(&((*scm).cmobile), &((*scs).cmobile));
  transcribe_array_int(&((*scm).cmobile_map), &((*scs).cmobile_map));
  transcribe_array_int(&((*scm).c_map), &((*scs).c_map));
  (*scs).solver = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj, (*scs).n_vars);
  (*scs).solver_data.f = Lagrange_solver_f;
  (*scs).solver_data.df = Lagrange_solver_df;
  (*scs).solver_data.fdf = Lagrange_solver_fdf;
  (*scs).solver_data.n = (*scs).n_vars;
  (*scs).solver_data.params = scs;
  (*scs).c_data = gsl_vector_alloc((*scs).n_vars);
  double *x0 = gsl_vector_ptr((*scs).c_data, 0);
  double *x0_ = gsl_vector_ptr((*scm).solver->x, 0);
  for (int mi = 0; mi < (*scm).n_vars; mi++)
    {
      x0[mi] = x0_[mi];
    }
  double L_init = total_length_func_exp((*scm).solver->x, (*scm).sc, &((*scm).c_map), &((*scm).top), &((*scm).cmobile), &((*scm).cmobile_map), 1e-32);
  heuristic_pars *hpars = (heuristic_pars *) (*scm).pars;
  x0[(*scs).tau_addr] = (*hpars).k * (L_init - L);
  double len = sc_solver_compute_fixed_length_exp(sc, &((*scs).top), &((*scs).cmobile), &((*scs).cmobile_map), &((*scs).c_map), (*scs).c_data);
  (*scs).len_disc = (*scs).L - len;
  gsl_multiroot_fdfsolver_set((*scs).solver, &((*scs).solver_data), (*scs).c_data);
}

void free_sc_Lagrange_solver(sc_Lagrange_solver *scs)
{
  free_contr_nbrlist(&((*scs).top));
  free_array_int(&((*scs).cmobile));
  free_array_int(&((*scs).cmobile_map));
  gsl_vector_free((*scs).c_data);
  gsl_multiroot_fdfsolver_free((*scs).solver);
  free_array_int(&((*scs).c_map));
}

void sc_Lagrange_solver_set(sc_Lagrange_solver *scs, int ext_f_site, double *ext_f, double L, double tol)
{
  if (ext_f_site > -1) (*scs).ext_f_site = ext_f_site;
  if (ext_f != NULL) for (int di = 0; di < (*(*scs).sc).dim; di++) (*scs).ext_f[di] = ext_f[di];
  if (L > 0) (*scs).L = L;
  (*scs).tol = tol;
  gsl_vector_memcpy((*scs).c_data, (*scs).solver->x);
  gsl_multiroot_fdfsolver_set((*scs).solver, &((*scs).solver_data), (*scs).c_data);
}

int sc_Lagrange_solver_relax_fdf_diag(sc_Lagrange_solver *scs, int N_steps, FILE *ofile)
{
  for (int ii = 0; ii < N_steps; ii++)
    {
      printf("Iteration %d:\n", ii);
      int status = gsl_multiroot_fdfsolver_iterate((*scs).solver);
      double *x = gsl_vector_ptr((*scs).solver-> x, 0);
      fprintf(ofile, "%d ", status);
      for (int i = 0; i < ((*scs).solver -> x) -> size; i++)
	{
	  fprintf(ofile, "%g ", x[i]);
	}
      gsl_vector *g = gsl_multiroot_fdfsolver_f((*scs).solver);
      fprintf(ofile, " %g\n", gsl_blas_dnrm2(g));
    }
}

int sc_Lagrange_solver_relax_fdf(sc_Lagrange_solver *scs, double tol)
{
  int count = 0;
  int status = 0;
  int status1, status2;
  while (count < MAX_ITER && status == 0)
    {
      count += 1;
      status = gsl_multiroot_fdfsolver_iterate((*scs).solver);
      if (status == GSL_ENOPROGJ)
	{
	  printf("Status = GSL_ENOPROGJ after %d iterations\n", count);
	  break;
	}
      gsl_vector *dx = gsl_multiroot_fdfsolver_dx((*scs).solver);
      double norm_dx = gsl_blas_dnrm2(dx);
      gsl_vector *f = gsl_multiroot_fdfsolver_f((*scs).solver);
      if (norm_dx < tol)
	{
	  // RESUME: consider including two tolerances in the Lagrange solver:
	  //           one for displacements, the other for the gradient
	  status1 = gsl_multiroot_test_residual(f, tol);
	  if (status1 == GSL_SUCCESS)
	    {
	      status = GSL_SUCCESS;
	      break;
	    }
	}
    }
  return status;
}

// NOTE: determine a lower bound on acceptable tolerances/precision
int string_config_relax_Lagrange(string_config *sc, int ext_f_site, double *ext_f, double L, double tol)
{
  sc_minimizer_mm scm;
  sc_minimizer_mm_init_heuristic(&scm, sc, ext_f_site, ext_f, L, 1.0);
  heuristic_pars *hpars = (heuristic_pars *) scm.pars;
  double stepsize = 0.002 * sqrt(string_config_var_x(sc));
  sc_minimizer_mm_solve(&scm, 1e-8);
  for (int i = 0; i < 2; i++)
    {
      gsl_vector_memcpy(scm.c_data, scm.solver->x);
      (*hpars).k *= 5;
      gsl_multimin_fdfminimizer_set(scm.solver, &(scm.solver_data), scm.c_data, stepsize, 1e-2);
      sc_minimizer_mm_solve(&scm, 1e-8);
      stepsize *= 0.2;
    }
  sc_Lagrange_solver scs;
  sc_Lagrange_solver_from_minimizer(&scs, &scm, ext_f_site, ext_f, L);
  free_sc_minimizer_mm(&scm);
  int status = sc_Lagrange_solver_relax_fdf(&scs, tol);
  if (status == GSL_SUCCESS)
    {
      printf("Transcribing coordinates\n");
      sc_solver_write_coords_exp(sc, scs.solver->x, &(scs.c_map));
    }
  else
    {
      printf("Lagrange relaxation halted with status code %d\n", status);
    }
  free_sc_Lagrange_solver(&scs);
}

// sc_cl_composite_solver

void sc_cl_composite_solver_init(sc_cl_composite_solver *sccs, string_config *sc, int ext_f_site, double *ext_f, double L, double tol)
{
  (*sccs).tol = tol;
  sc_minimizer_mm_init_heuristic(&((*sccs).scm), sc, ext_f_site, ext_f, L, 1.0);
  (*sccs).hpars = (heuristic_pars *) (*sccs).scm.pars;
  int status1 = sc_minimizer_mm_solve(&((*sccs).scm), tol);
  sc_Lagrange_solver_from_minimizer(&((*sccs).scLs), &((*sccs).scm), ext_f_site, ext_f, L);
  int status2 = sc_Lagrange_solver_relax_fdf(&((*sccs).scLs), tol);
}

void free_sc_cl_composite_solver(sc_cl_composite_solver *sccs)
{
  free_sc_minimizer_mm(&((*sccs).scm));
  free_sc_Lagrange_solver(&((*sccs).scLs));
}

double sc_cl_composite_solver_iterate(sc_cl_composite_solver *sccs)
{
  transcribe_contr_nbrlist(&((*sccs).scm.top), &((*sccs).ref_top));
  transcribe_array_int(&((*sccs).scm.cmobile), &((*sccs).ref_cmobile));
  transcribe_array_int(&((*sccs).scm.cmobile_map), &((*sccs).ref_cmobile_map));
  transcribe_array_int(&((*sccs).scm.c_map), &((*sccs).ref_c_map));
  (*sccs).ref_data = gsl_vector_alloc((*sccs).scLs.n_vars);
  gsl_vector_memcpy((*sccs).ref_data, (*sccs).scLs.solver->x);
  ((*sccs).hpars->k) *= 2;
  sc_minimizer_mm_solve(&((*sccs).scm), (*sccs).tol);
  free_sc_Lagrange_solver(&((*sccs).scLs));
  sc_Lagrange_solver_from_minimizer(&((*sccs).scLs), &((*sccs).scm), (*sccs).hpars->ext_f_site, (*sccs).hpars->ext_f, (*sccs).hpars->L0);
  int status = sc_Lagrange_solver_relax_fdf(&((*sccs).scLs), (*sccs).tol);
  double distsq = sc_solver_distsq_exp((*sccs).scLs.sc, &((*sccs).ref_top), &((*sccs).ref_cmobile), &((*sccs).ref_cmobile_map), &((*sccs).ref_c_map), (*sccs).ref_data, &((*sccs).scLs.top), &((*sccs).scLs.cmobile), &((*sccs).scLs.cmobile_map), &((*sccs).scLs.c_map), (*sccs).scLs.solver->x);
  free_contr_nbrlist(&((*sccs).ref_top));
  free_array_int(&((*sccs).ref_cmobile));
  free_array_int(&((*sccs).ref_cmobile_map));
  free_array_int(&((*sccs).ref_c_map));
  gsl_vector_free((*sccs).ref_data);
  return sqrt(distsq);
}					 

// RESUME: test this!
void sc_solver_expand(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data, gsl_vector *nc_data)
{
  int init_len_top = (*top).top.top.v.len;
  int ci = init_len_top;
  int Nm0 = (*cmobile).len;
  while (ci > 0)
    {
      ci -= 1;
      if ((*top).fibers.e[ci].len > 1)
	{
	  for (int fi = 1; fi < (*top).fibers.e[ci].len; fi++)
	    {
	      add2array_int(c_map, (*top).fibers.e[ci].e[0]);
	    }
	  sc_solver_expand_cluster_exp(sc, top, cmobile, cmobile_map, ci);
	}
    }
  if ((*sc).mobile.len == (*cmobile).len) {}
  else
    {
      printf("Something weird happened! a;lsdf;oajwa;ea;\n");
      exit(EXIT_FAILURE);
    }
  int base = 0;
  for (int mi = 0; mi < Nm0; mi++)
    {
      double *xi = gsl_vector_ptr(nc_data, base);
      double *xi_ = gsl_vector_ptr(c_data, base);
      for (int di = 0; di < (*sc).dim; di++) xi[di] = xi_[di];
      base += (*sc).dim;
    }
  for (int mi = Nm0; mi < (*cmobile).len; mi++)
    {
      double *xi = gsl_vector_ptr(nc_data, base);
      double *xi_;
      int ri = (*c_map).e[mi];
      int ci = (*top).map.e[ri];
      if (ci < init_len_top) {}
      else
	{
	  printf("Something quite strange occurred... BLEARGGH\n");
	  exit(EXIT_FAILURE);
	}
      if ((*cmobile_map).e[ci] > -1)
	{
	  // Check that ci is consistent with assumptions: recall 'fi' is fibers.e[ci].e[0] for some ci < init_len(top)
	  int mci = (*cmobile_map).e[ci];
	  xi_ = gsl_vector_ptr(c_data, (*c_map).e[mci]);
	}
      else
	{
	  xi_ = string_config_vertex_coords(sc, ri);
	}
      for (int di = 0; di < (*sc).dim; di++) xi[di] = xi_[di];
      (*c_map).e[mi] = base;
      base += (*sc).dim;
    }
}

// Graph operations

void contract_leaves(contr_nbrlist *aux, array_int *bdry, array_int *mobile, array_int *is_mobile)
{
  if ((*mobile).len > 0) {}
  else return;
  array_int tmp;
  tmp.mem = 0;
  if (bdry != NULL) {}
  else
    {
      bdry = &tmp;
      array_int_init(&tmp, 0);
      for (int mi = 0; mi < (*mobile).len; mi++)
	{
	  int ci = (*mobile).e[mi];
	  if ((*is_mobile).e[ci] > -1) {}
	  else
	    {
	      printf("Inconsistency found in 'mobility' array/map\n");
	      exit(EXIT_FAILURE);
	    }
	  if ((*aux).top.top.v.e[ci].len == 1) add2array_int(&tmp, ci);
	}
    }
  if ((*bdry).len > 0)
    {
      nbrlist *top = &((*aux).top.top);
      // RESUME: implement this, change contract_elbows for new 'mobile' structure
      array_int bdry_addr; // NOTE: consider replacing this with a hash table
      array_int_init(&bdry_addr, (*top).v.len);
      bdry_addr.len = (*top).v.len;
      for (int i = 0; i < (*top).v.len; i++) bdry_addr.e[i] = -1;
      for (int i = 0; i < (*bdry).len; i++) bdry_addr.e[(*bdry).e[i]] = i;
      while ((*bdry).len > 0)
	{
	  int blenm1 = (*bdry).len - 1;
	  int ci = (*bdry).e[blenm1];
	  remove_array_int(bdry, blenm1);
	  if ((*aux).top.top.v.e[ci].len == 1 && (*is_mobile).e[ci] > -1) // This should hold trivially
	    { 
	      int cii = (*aux).top.top.v.e[ci].e[0];
	      if ((*is_mobile).e[cii] > -1 && (*aux).top.top.v.e[cii].len == 2)
		{
		  bdry_addr.e[cii] = (*bdry).len;
		  add2array_int(bdry, cii);
		}
	      printf("Attempting to contract %d to %d\n", ci, cii);
	      contr_nbrlist_merge_cluster(aux, ci, cii);
	      remove_array_int(&bdry_addr, ci); // This should map the last element to 'ci'
	      if (ci < bdry_addr.len && bdry_addr.e[ci] > -1) (*bdry).e[bdry_addr.e[ci]] = ci;
	      int mi = (*is_mobile).e[ci];
	      remove_array_int(is_mobile, ci); // This should map the last element to 'ci'
	      if (ci < (*is_mobile).len && (*is_mobile).e[ci] > -1) (*mobile).e[(*is_mobile).e[ci]] = ci;
	      remove_array_int(mobile, mi);
	      if (mi < (*mobile).len)
		{
		  (*is_mobile).e[(*mobile).e[mi]] = mi;
		  if ((*mobile).e[mi] < (*aux).top.top.v.len) {}
		  else
		    {
		      printf("Error: inconsistency between mobile array and contracted topology (%d %d %d)\n", (*mobile).e[mi], (*aux).top.top.v.len, (*is_mobile).len);
		      exit(EXIT_FAILURE);
		    }
		}
	    }
	  else
	    {
	      if (bdry_addr.e[ci] == blenm1) {}
	      else printf("(Inconsistency in bdry list\n");
	      printf("Something weird happened! oawieur892 n_nbors = %d, mobile_addr = %d\n", (*aux).top.top.v.e[ci].len, (*is_mobile).e[ci]);
	      exit(EXIT_FAILURE);
	    }
	}
      free_array_int(&bdry_addr);
      if (tmp.mem > 0)
	{
	  free_array_int(&tmp);
	}
    }
}

void contract_elbows(contr_nbrlist *aux, array_int *elbows, array_int *mobile, array_int *is_mobile)
{
  if ((*mobile).len > 0) {}
  else return;
  //printf("contract_elbows\n");
  array_int tmp;
  tmp.mem = 0;
  if (elbows == NULL)
    {
      //printf("\tFinding elbows\n");
      elbows = &tmp;
      array_int_init(&tmp, 0);
      //      if (mobile == NULL) printf("Warning: mobile = null\n");
      for (int mi = 0; mi < (*mobile).len; mi++)
	{
	  int ci = (*mobile).e[mi];
	  if ((*aux).top.top.v.e[ci].len == 2) add2array_int(&tmp, ci);
	}
      //printf("\t(done)\n");
    }
  if ((*elbows).len > 0)
    {
      //printf("\tNontrivial elbows found\n");
      array_int elbows_addr;
      array_int_init(&elbows_addr, (*aux).top.top.v.len);
      elbows_addr.len = (*aux).top.top.v.len;
      for (int i = 0; i < (*aux).top.top.v.len; i++) elbows_addr.e[i] = -1;
      for (int ei = 0; ei < (*elbows).len; ei++) elbows_addr.e[(*elbows).e[ei]] = ei;
      while ((*elbows).len > 0)
	{
	  int ci = (*elbows).e[0];
	  if ((*aux).top.top.v.e[ci].len == 2)
	    {
	      int cii = (*aux).top.top.v.e[ci].e[0];
	      int ciii = (*aux).top.top.v.e[ci].e[1];
	      int wii = (*aux).top.edge_wts.e[ci].e[0];
	      int wiii = (*aux).top.edge_wts.e[ci].e[1];
	      if (wii >= wiii) contr_nbrlist_merge_cluster(aux, ci, cii);
	      else contr_nbrlist_merge_cluster(aux, ci, ciii);
	      int mi = (*is_mobile).e[ci];
	      remove_array_int(mobile, mi);
	      remove_array_int(is_mobile, ci);
	      if ((*mobile).len > mi) (*is_mobile).e[(*mobile).e[mi]] = mi;
	      if ((*is_mobile).len > ci && (*is_mobile).e[ci] > -1) (*mobile).e[(*is_mobile).e[ci]] = ci;
	    }
	  remove_array_int(elbows, 0);
	  if ((*elbows).len > 0) elbows_addr.e[(*elbows).e[0]] = 0;
	  remove_array_int(&elbows_addr, ci);
	  if (elbows_addr.len > ci && elbows_addr.e[ci] > -1) (*elbows).e[elbows_addr.e[ci]] = ci;
	}
      free_array_int(&elbows_addr);
      if (tmp.mem > 0)
	{
	  free_array_int(&tmp);
	}
     }
 }

// 'prep' a contractable list (of integers)
void prep_contr_list(array_int *src, array_int *dest, array_int *dest_map, int map_size)
{
  array_int_init(dest, (*src).len);
  array_int_init(dest_map, map_size);
  (*dest_map).len = map_size;
  (*dest).len = (*src).len;
  for (int i = 0; i < map_size; i++) (*dest_map).e[i] = -1;
  for (int mi = 0; mi < (*src).len; mi++)
    {
      int i = (*src).e[mi];
      (*dest).e[mi] = i;
      (*dest_map).e[i] = mi;
    }
}

double sc_minimizer_shortest_dist(string_config *sc, contr_nbrlist *top, array_int *cmobile, array_int *cmobile_map, array_int *c_map, gsl_vector *c_data)
{
  double min_distsq = 1e99;
  //if ((*sc).mobile.len > (*cmobile).len) return 0;
  for (int ci = 0; ci < (*top).top.top.v.len; ci++)
    {
      int mci = (*cmobile_map).e[ci];
      double *xci;
      if (mci > -1) xci = gsl_vector_ptr(c_data, (*c_map).e[mci]);
      else
	{
	  int rci = (*top).fibers.e[ci].e[0];
	  xci = string_config_vertex_coords(sc, rci);
	}
      for (int ni = 0; ni < (*top).top.top.v.e[ci].len; ni++)
	{
	  int cii = (*top).top.top.v.e[ci].e[ni];
	  if (cii < ci) {}
	  else continue;
	  int mcii = (*cmobile_map).e[cii];
	  double *xcii;
	  if (mcii > -1)
	    {
	      xcii = gsl_vector_ptr(c_data, (*c_map).e[mcii]);
	    }
	  else
	    {
	      int rcii = (*top).fibers.e[cii].e[0];
	      xcii = string_config_vertex_coords(sc, rcii);
	    }
	  double distsq = euclid_distsq(xci, xcii, (*sc).dim);
	  if (distsq < min_distsq) min_distsq = distsq;
	  if (min_distsq == 0) break;
	}
    }
  return sqrt(min_distsq);
}

// Arrangements with multiple strings
// Tack sets

void tack_set_init(tack_set *ts, int dim)
{
  (*ts).dim = dim;
  array_voidstar_init(&((*ts).coords), 0);
  array_char_init(&((*ts).is_fxd), 0);
  array_int_init(&((*ts).mobile), 0);
  array_int_init(&((*ts).fxd), 0);
  array_int_init(&((*ts).fm_addr), 0);
}

double tack_set_var_x(tack_set *ts)
{
  double ave[(*ts).dim];
  for (int di = 0; di < (*ts).dim; di++) ave[di] = 0;
  for (int i = 0; i < (*ts).coords.len; i++)
    {
      double *xi = (double *) (*ts).coords.e[i];
      for (int di = 0; di < (*ts).dim; di++) ave[di] += xi[di];
    }
  double wt = 1. / (*ts).coords.len;
  for (int di = 0; di < (*ts).dim; di++) ave[di] *= wt;
  double var = 0;
  for (int i = 0; i < (*ts).coords.len; i++)
    {
      double *xi = (double *) (*ts).coords.e[i];
      for (int di = 0; di < (*ts).dim; di++)
	{
	  double delx = xi[di] - ave[di];
	  var += delx * delx;
	}
    }
  return var / (*ts).coords.len;
}

void add_fxd_pt_tack_set(tack_set *ts, double *x)
{
  add2array_int(&((*ts).fm_addr), (*ts).fxd.len);
  add2array_int(&((*ts).fxd), (*ts).is_fxd.len);
  add2array_char(&((*ts).is_fxd), 1);
  add2array_voidstar(&((*ts).coords), x);
}

void add_mbl_pt_tack_set(tack_set *ts, double *x)
{
  add2array_int(&((*ts).fm_addr), (*ts).mobile.len);
  add2array_int(&((*ts).mobile), (*ts).is_fxd.len);
  add2array_char(&((*ts).is_fxd), 0);
  add2array_voidstar(&((*ts).coords), x);
}

void free_tack_set(tack_set *ts, void (*free_coord_elem)(void *))
{
  free_array_int(&((*ts).fxd));
  free_array_int(&((*ts).mobile));
  free_array_int(&((*ts).fm_addr));
  free_array_char(&((*ts).is_fxd));
  free_array_voidstar(&((*ts).coords), free_coord_elem);
}

void free_tack_set_shallow(tack_set *ts)
{
  free_array_int(&((*ts).fxd));
  free_array_int(&((*ts).mobile));
  free_array_int(&((*ts).fm_addr));
  free_array_char(&((*ts).is_fxd));
  free_array_voidstar(&((*ts).coords), NULL);
}

void linked_sc_init(linked_sc *lsc, tack_set *ts)
{
  (*lsc).ts = ts;
  array_voidstar_init(&((*lsc).tops), 0);
  array_voidstar_init(&((*lsc).embs), 0);
}

char linked_sc_connected(linked_sc *lsc)
{
  array_int uf;
  union_find_bb_init(&uf, (*lsc).ts->coords.len);
  for (int ti = 0; ti < (*lsc).tops.len; ti++)
    {
      edge_wtd_graph *top = (edge_wtd_graph *) (*lsc).tops.e[ti];
      array_int *emb = (array_int *) (*lsc).embs.e[ti];
      for (int tti = 0; tti < (*top).top.v.len; tti++)
	{
	  int itti = (*emb).e[tti];
	  for (int ni = 0; ni < (*top).top.v.e[tti].len; ni++)
	    {
	      int ttii = (*top).top.v.e[tti].e[ni];
	      int ittii = (*emb).e[ttii];
	      union_find_bb_union(&uf, itti, ittii);
	    }
	}
    }
  int r0;
  union_find_bb_find(&uf, 0, &r0);
  int status = 1;
  for (int i = 1; i < (*lsc).ts->coords.len; i++)
    {
      int ri;
      union_find_bb_find(&uf, i, &ri);
      if (ri == r0) {}
      else
	{
	  printf("linked_sc is disconnected\n");
	  status = 0;
	  break;
	}
    }
  free_array_int(&uf);
  return status;
}

void free_linked_sc_shallow(linked_sc *lsc)
{
  free_array_voidstar(&((*lsc).tops), NULL);
  free_array_voidstar(&((*lsc).embs), NULL);
}

void free_linked_sc(linked_sc *lsc)
{
  for (int i = 0; i < (*lsc).tops.len; i++)
    {
      edge_wtd_graph *top = (edge_wtd_graph *) (*lsc).tops.e[i];
      array_int *emb = (array_int *) (*lsc).embs.e[i];
      free_edge_wtd_graph(top);
      free_array_int(emb);
    }
  free_linked_sc_shallow(lsc);
}

void add2linked_sc(linked_sc *lsc, edge_wtd_graph *top, array_int *emb)
{
  add2array_voidstar(&((*lsc).tops), top);
  add2array_voidstar(&((*lsc).embs), emb);
}

double contr_pointset_distsq(const gsl_vector *x0, coord_pars *cpars0, const gsl_vector *x1, coord_pars *cpars1)
{
 tack_set *ts = (*cpars0).ts;
  double distsq = 0;
  for (int mi = 0; mi < (*ts).mobile.len; mi++)
    {
      int vi = (*ts).mobile.e[mi];
      int cvi0 = (*cpars0).map.e[vi];
      int cvi1 = (*cpars1).map.e[vi];
      int mcvi0 = (*cpars0).cmobile_map.e[cvi0];
      int mcvi1 = (*cpars1).cmobile_map.e[cvi1];
      const double *xvi0;
      const double *xvi1;
      if (mcvi0 > -1) xvi0 = gsl_vector_const_ptr(x0, (*cpars0).c_map.e[mcvi0]);
      else xvi0 = (double *) (*ts).coords.e[(*cpars0).fibers.e[cvi0].e[0]];
      if (mcvi1 > -1) xvi1 = gsl_vector_const_ptr(x1, (*cpars1).c_map.e[mcvi1]);
      else xvi1 = (double *) (*ts).coords.e[(*cpars1).fibers.e[cvi1].e[0]];
      distsq += euclid_distsq(xvi0, xvi1, (*ts).dim);
    }
  return distsq;
}

double contr_pointset_dist(const gsl_vector *x0, coord_pars *cpars0, const gsl_vector *x1, coord_pars *cpars1)
{
  return sqrt(contr_pointset_distsq(x0, cpars0, x1, cpars1));
}

void coord_pars_from_tack_set(coord_pars *cpars, tack_set *ts)
{
  (*cpars).ts = ts;
  prep_contr_list(&((*ts).mobile), &((*cpars).cmobile), &((*cpars).cmobile_map), (*ts).coords.len);
  array_int_init(&((*cpars).map), (*ts).coords.len);
  aarray_int_init_precise(&((*cpars).fibers), (*ts).coords.len, 1);
  array_int_init(&((*cpars).c_map), (*ts).mobile.len);
  for (int i = 0; i < (*ts).coords.len; i++)
    {
      extend_aarray_int(&((*cpars).fibers));
      add2array_int(&((*cpars).fibers.e[i]), i);
      add2array_int(&((*cpars).map), i);
    }
  int base = 0;
  for (int mi = 0; mi < (*ts).mobile.len; mi++)
    {
      (*cpars).c_map.e[mi] = base;
      base += (*ts).dim;
    }
}

void transcribe_coord_pars(coord_pars *src, coord_pars *dest)
{
  (*dest).ts = (*src).ts;
  transcribe_array_int(&((*src).cmobile), &((*dest).cmobile));
  transcribe_array_int(&((*src).cmobile_map), &((*dest).cmobile_map));
  transcribe_array_int(&((*src).c_map), &((*dest).c_map));
  transcribe_aarray_int(&((*src).fibers), &((*dest).fibers));
  transcribe_array_int(&((*src).map), &((*dest).map));
}

void free_coord_pars(coord_pars *cpars)
{
  free_array_int(&((*cpars).cmobile));
  free_array_int(&((*cpars).cmobile_map));
  free_array_int(&((*cpars).c_map));
  free_aarray_int(&((*cpars).fibers));
  free_array_int(&((*cpars).map));
}

// Merge cluster ci into cluster cj without updating c_data
void coord_pars_merge(coord_pars *cpars, int ci, int cj)
{
  // fibers, map, cmobile, cmobile_map, c_map
  for (int fi = 0; fi < (*cpars).fibers.e[ci].len; fi++)
    {
      int i = (*cpars).fibers.e[ci].e[fi];
      add2array_int(&((*cpars).fibers.e[cj]), i);
      (*cpars).map.e[i] = cj;
    }
  remove_aarray_int(&((*cpars).fibers), ci);
  if (ci < (*cpars).fibers.len)
    {
      for (int fi = 0; fi < (*cpars).fibers.e[ci].len; fi++)
	{
	  int i = (*cpars).fibers.e[ci].e[fi];
	  (*cpars).map.e[i] = ci;
	}
    }
  int mci = (*cpars).cmobile_map.e[ci];
  if (mci > -1)
    {
      remove_array_int(&((*cpars).cmobile), mci);
      remove_array_int(&((*cpars).c_map), mci);
      if (mci < (*cpars).cmobile.len)
	{
	  int cii = (*cpars).cmobile.e[mci];
	  (*cpars).cmobile_map.e[cii] = mci;
	}
    }
  remove_array_int(&((*cpars).cmobile_map), ci);
  if (ci < (*cpars).cmobile_map.len)
    {  
      int mci_ = (*cpars).cmobile_map.e[ci];
      if (mci_ > -1) (*cpars).cmobile.e[mci_] = ci;
    }
}

void lsc_heuristic_pars_init(lsc_heuristic_pars *hpars, linked_sc *lsc, int i, double *ext_f, double *Ls, double k, coord_pars *cpars, array_voidstar *ctops, array_voidstar *cembs, aarray_int *cembs_map)
{
  (*hpars).lsc = lsc;
  (*hpars).i = i;
  (*hpars).ext_f = ext_f;
  (*hpars).Ls = Ls;
  (*hpars).k = k;
  (*hpars).cpars = cpars;
  (*hpars).ctops = ctops;
  (*hpars).cembs = cembs;
  (*hpars).cembs_map = cembs_map;
}

double lsc_solver_var_coords(coord_pars *cpars, const gsl_vector *x)
{
  gsl_vector *ave_coords = gsl_vector_alloc((*(*cpars).ts).dim);
  gsl_vector_set_zero(ave_coords);
  double *ax = gsl_vector_ptr(ave_coords, 0);
  for (int ci = 0; ci < (*cpars).fibers.len; ci++)
    {
      int mci = (*cpars).cmobile_map.e[ci];
      const double *xci;
      if (mci > -1) xci = gsl_vector_const_ptr(x, (*cpars).c_map.e[mci]);
      else
	{
	  int rci = (*cpars).fibers.e[ci].e[0];
	  xci = (double *) (*(*cpars).ts).coords.e[rci];
	}
      for (int di = 0; di < (*(*cpars).ts).dim; di++) ax[di] += xci[di];
    }
  double wt = 1. / (*cpars).fibers.len;
  gsl_vector_scale(ave_coords, wt);
  double var_coords = 0;
  for (int ci = 0; ci < (*cpars).fibers.len; ci++)
    {
      int mci = (*cpars).cmobile_map.e[ci];
      const double *xci;
      if (mci > -1) xci = gsl_vector_const_ptr(x, (*cpars).c_map.e[mci]);
      else
	{
	  int rci = (*cpars).fibers.e[ci].e[0];
	  xci = (double *) (*(*cpars).ts).coords.e[rci];
	}
      var_coords += euclid_distsq(xci, &ax[0], (*(*cpars).ts).dim);
    }
  var_coords *= wt;
  return var_coords;
}

void lsc_solver_transcribe_tops(array_voidstar *ctops, array_voidstar *tops)
{
  array_voidstar_init(ctops, (*tops).len);
  for (int i = 0; i < (*tops).len; i++)
    {
      edge_wtd_graph *ewg = (edge_wtd_graph *) malloc(sizeof(edge_wtd_graph));
      transcribe_edge_wtd_graph((edge_wtd_graph *) (*tops).e[i], ewg);
      add2array_voidstar(ctops, ewg);
    }
}

void lsc_solver_transcribe_embs(array_voidstar *cembs, aarray_int *cembs_map, array_voidstar *embs, int range_size)
{
  array_voidstar_init(cembs, (*embs).len);
  aarray_int_init_precise(cembs_map, range_size, (*embs).len);
  (*cembs_map).len = range_size;
  for (int i = 0; i < range_size; i++)
    {
      (*cembs_map).e[i].len = (*embs).len;
      for (int ti = 0; ti < (*embs).len; ti++) (*cembs_map).e[i].e[ti] = -1;
    }
  for (int ti = 0; ti < (*embs).len; ti++)
    {
      array_int *emb = (array_int *) malloc(sizeof(array_int));
      transcribe_array_int((array_int *) (*embs).e[ti], emb);
      add2array_voidstar(cembs, emb);
      for (int ei = 0; ei < (*emb).len; ei++) (*cembs_map).e[(*emb).e[ei]].e[ti] = ei;
    }
}

void lsc_h_minimizer_init(lsc_h_minimizer *lscm, linked_sc *lsc, int i, double *ext_f, double *Ls, double k)
{
  if (ext_f != NULL)
    {
      if ((*lsc).ts->is_fxd.e[i])
	{
	  printf("Error: attempting to address immobile vertex\n");
	  exit(EXIT_FAILURE);
	}
    }
  (*lscm).i = i;
  (*lscm).ext_f = ext_f;
  tack_set *ts = (*lsc).ts;
  (*lscm).Ls = Ls;
  (*lscm).lsc = lsc;
  lsc_solver_transcribe_tops(&((*lscm).ctops), &((*lsc).tops));
  lsc_solver_transcribe_embs(&((*lscm).cembs), &((*lscm).cembs_map), &((*lsc).embs), (*ts).coords.len);
  coord_pars_from_tack_set(&((*lscm).cpars), (*lsc).ts);
  lsc_heuristic_pars_init(&((*lscm).hpars), lsc,  i, ext_f, Ls, k, &((*lscm).cpars), &((*lscm).ctops), &((*lscm).cembs), &((*lscm).cembs_map));
  (*lscm).hsolver_data.f = lsc_h_f;
  (*lscm).hsolver_data.df = lsc_h_df;
  (*lscm).hsolver_data.fdf = lsc_h_fdf;
  (*lscm).hsolver_data.params = &((*lscm).hpars);
  (*lscm).hsolver_data.n = (*(*lsc).ts).mobile.len * (*(*lsc).ts).dim;
  (*lscm).hsolver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, (*lscm).hsolver_data.n);
  (*lscm).c_data = gsl_vector_alloc((*lscm).hsolver_data.n);
  // Load mobile coordinates into c_data
  //printf("Defining coordinates:\n");
  for (int mi = 0; mi < (*ts).mobile.len; mi++)
    {
      double *xi = gsl_vector_ptr((*lscm).c_data, (*lscm).cpars.c_map.e[mi]);
      double *xi_ = (double *) (*ts).coords.e[(*ts).mobile.e[mi]];
      for (int di = 0; di < (*ts).dim; di++) xi[di] = xi_[di];
    }
  double *c_data_ = gsl_vector_ptr((*lscm).c_data, 0);
  double stepsize = sqrt(tack_set_var_x(ts));
  for (int i = 0; i < (*(*lscm).c_data).size; i++) c_data_[i] += stepsize * (rnd() - 0.5);
  gsl_multimin_fdfminimizer_set((*lscm).hsolver, &((*lscm).hsolver_data), (*lscm).c_data, stepsize, 1e-2);
  for (int i = 0; i < 3; i++) lsc_h_minimizer_iterate(lscm);
}

void free_lsc_h_minimizer(lsc_h_minimizer *lscm)
{
  gsl_vector_free((*lscm).c_data);
  gsl_multimin_fdfminimizer_free((*lscm).hsolver);
  free_array_voidstar(&((*lscm).ctops), free_edge_wtd_graph);
  free_array_voidstar(&((*lscm).cembs), free_array_int);
  free_aarray_int(&((*lscm).cembs_map));
  free_coord_pars(&((*lscm).cpars));
}

int lsc_h_minimizer_iterate(lsc_h_minimizer *lscm)
{
  return gsl_multimin_fdfminimizer_iterate((*lscm).hsolver);
}

// RESUME: check this!
void lsc_h_minimizer_check_merging(lsc_h_minimizer *lscm, double merge_rad)
{
  int init_n_mobile = (*lscm).cpars.cmobile.len;
  double mrsq = merge_rad * merge_rad;
  int mci = 0;
  while (mci < (*lscm).cpars.cmobile.len)
    {
      int ci = (*lscm).cpars.cmobile.e[mci];
      double *xci = gsl_vector_ptr((*lscm).hsolver->x, (*lscm).cpars.c_map.e[mci]);
      int cii_ = -1;
      for (int ti = 0; ti < (*lscm).ctops.len; ti++)
	{
	  if ((*lscm).cembs_map.e[ci].e[ti] > -1)
	    {
	      array_int *cemb_ti = (array_int *) (*lscm).cembs.e[ti];
	      edge_wtd_graph *ewg = (edge_wtd_graph *) (*lscm).ctops.e[ti];
	      int iti = (*lscm).cembs_map.e[ci].e[ti];
	      for (int ni = 0; ni < (*ewg).top.v.e[iti].len; ni++)
		{
		  int iiti = (*ewg).top.v.e[iti].e[ni];
		  int cii = (*cemb_ti).e[iiti];
		  int mcii = (*lscm).cpars.cmobile_map.e[cii];
		  double *xcii;
		  if (mcii > -1) xcii = gsl_vector_ptr((*lscm).hsolver->x, (*lscm).cpars.c_map.e[mcii]);
		  else
		    {
		      int rii = (*lscm).cpars.fibers.e[cii].e[0];
		      xcii = (double *) ((*lscm).lsc->ts->coords).e[rii];
		    }
		  double distsq = euclid_distsq(xci, xcii, (*lscm).lsc->ts->dim);
		  if (distsq < mrsq)
		    {
		      cii_ = cii;
		      break;
		    }
		}
	      if (cii_ > -1) break;
	    } 
	}
      if (cii_ > -1)
	{
	  // Merge the associated clusters in ctops, and adjust the associated embeddings
	  for (int tii = 0; tii < (*lscm).ctops.len; tii++)
	    {
	      edge_wtd_graph *top = (edge_wtd_graph *) (*lscm).ctops.e[tii];
	      array_int *emb = (array_int *) (*lscm).cembs.e[tii];
	      int tiici = (*lscm).cembs_map.e[ci].e[tii];
	      int tiicii_ = (*lscm).cembs_map.e[cii_].e[tii];
	      if (tiici > -1 && tiicii_ > -1)
		{
		  edge_wtd_graph_merge(top, tiici, tiicii_);
		  remove_array_int(emb, tiici);
		  if ((*emb).len > tiici)
		    {
		      int ci_ = (*emb).e[tiici];
		      (*lscm).cembs_map.e[ci_].e[tii] = tiici;
		    }
		  (*lscm).cembs_map.e[ci].e[tii] = tiicii_;
		}
	      else if (tiici > -1)
		{
		  (*lscm).cembs_map.e[cii_].e[tii] = tiici;
		  (*emb).e[tiici] = cii_;
		}
	      else if (tiicii_ > -1)
		{
		  (*lscm).cembs_map.e[ci].e[tii] = tiicii_;
		}
	    }
	  coord_pars_merge(&((*lscm).cpars), ci, cii_);
	  // Update cembs_map and possibly cembs
	  remove_aarray_int(&((*lscm).cembs_map), ci);
	  if (ci < (*lscm).cembs_map.len)
	    {
	      for (int tii = 0; tii < (*lscm).ctops.len; tii++)
		{
		  array_int *emb_tii = (array_int *) (*lscm).cembs.e[tii];
		  int tiici = (*lscm).cembs_map.e[ci].e[tii];
		  if (tiici > -1) (*emb_tii).e[tiici] = ci;
		}
	    }
	}
      else mci += 1;
    }
  // Redefine coordinates (and c_map), reset the solver
  if (init_n_mobile == (*lscm).cpars.cmobile.len) return;
  int n_vars = (*lscm).cpars.cmobile.len * (*(*lscm).lsc).ts->dim;
  gsl_vector *nc_data = gsl_vector_alloc(n_vars);
  int base = 0;
  for (int mi = 0; mi < (*lscm).cpars.cmobile.len; mi++)
    {
      double *xi = gsl_vector_ptr(nc_data, base);
      double *xi_ = gsl_vector_ptr((*lscm).hsolver->x, (*lscm).cpars.c_map.e[mi]);
      for (int di = 0; di < (*(*lscm).lsc).ts->dim; di++) xi[di] = xi_[di];
      (*lscm).cpars.c_map.e[mi] = base;
      base += (*(*lscm).lsc).ts->dim;
    }
  gsl_vector_free((*lscm).c_data);
  (*lscm).c_data = nc_data;
  gsl_multimin_fdfminimizer_free((*lscm).hsolver);
  (*lscm).hsolver_data.n = n_vars;
  (*lscm).hsolver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
  double stepsize = 0.01 * sqrt(tack_set_var_x((*lscm).lsc->ts));
  gsl_multimin_fdfminimizer_set((*lscm).hsolver, &((*lscm).hsolver_data), (*lscm).c_data, stepsize, 1e-2);
}

void lsc_h_minimizer_expand_clusters(lsc_h_minimizer *lscm)
{
  // Reset the coord_pars and (contracted) topologies associated with lscm, expand coordinates to include original mobile vertices.
  tack_set *ts = (*lscm).lsc->ts;
  array_int *ts_mobile = &(ts->mobile);
  int n_vars = (*ts_mobile).len * ts->dim;
  gsl_vector *ncdata = gsl_vector_alloc(n_vars);
  int base = 0;
  for (int mi = 0; mi < (*ts_mobile).len; mi++)
    {
      int i = (*ts_mobile).e[mi];
      double *xi = gsl_vector_ptr(ncdata, base);
      double *xi_;
      int ci = (*lscm).cpars.map.e[i];
      int mci = (*lscm).cpars.cmobile_map.e[ci];
      if (mci > -1)
	{
	  xi_ = gsl_vector_ptr((*lscm).hsolver->x, (*lscm).cpars.c_map.e[mci]);
	}
      else
	{
	  int rci = (*lscm).cpars.fibers.e[ci].e[0];
	  xi_ = (double *) (*ts).coords.e[rci];
	}
      for (int di = 0; di < (*ts).dim; di++) xi[di] = xi_[di];
      base += (*ts).dim;
    }
  free_coord_pars(&((*lscm).cpars));
  coord_pars_from_tack_set(&((*lscm).cpars), ts); // This should also reset c_map
  free_array_voidstar(&((*lscm).ctops), free_edge_wtd_graph);
  free_array_voidstar(&((*lscm).cembs), free_array_int);
  free_aarray_int(&((*lscm).cembs_map));
  // Reset contracted topologies and embeddings from lsc
  linked_sc *lsc = (*lscm).lsc;
  lsc_solver_transcribe_tops(&((*lscm).ctops), &((*lsc).tops));
  lsc_solver_transcribe_embs(&((*lscm).cembs), &((*lscm).cembs_map), &((*lsc).embs), (*ts).coords.len);
  gsl_vector_free((*lscm).c_data);
  (*lscm).c_data = ncdata;
  (*lscm).hsolver_data.n = (*ncdata).size;
  double stepsize = sqrt(lsc_solver_var_coords(&((*lscm).cpars), ncdata)) * 1e-2;
  gsl_multimin_fdfminimizer_free((*lscm).hsolver);
  (*lscm).hsolver = gsl_multimin_fdfminimizer_alloc(gsl_minimizer_type_mm, n_vars);
  gsl_multimin_fdfminimizer_set((*lscm).hsolver, &((*lscm).hsolver_data), (*lscm).c_data, stepsize, 1e-2);
}

int lsc_h_minimizer_relax(lsc_h_minimizer *lscm, double tol)
{
  //  printf("lsc_h_minimizer_relax: %g\n", tol);
  // Iterate, choose a convergence criterion (e.g. gradient norm)
  double tolsq = tol * tol;
  double merge_rad = 5 * tol;
  int count = 0;
  int status = GSL_CONTINUE;
  while (count < MAX_ITER)
    {
      count += 1;
      for (int j = 0; j < 5; j++) status = lsc_h_minimizer_iterate(lscm);
      if (status == GSL_SUCCESS)
	{
	  gsl_vector *df = gsl_multimin_fdfminimizer_gradient((*lscm).hsolver);
	  gsl_vector *dx = gsl_multimin_fdfminimizer_dx((*lscm).hsolver);
	  double *df_ = gsl_vector_ptr(df, 0);
	  double *dx_ = gsl_vector_ptr(dx, 0);
	  double err = euclid_normsq(df_, (*(*lscm).lsc).ts->dim) + euclid_normsq(dx_, (*(*lscm).lsc).ts->dim);
	  if (err < tolsq) break;
	  else status = GSL_CONTINUE;
	}
      else break;
      lsc_h_minimizer_check_merging(lscm, merge_rad);
    }
  //printf("(done) count = %d, status = %d, dof = %ld\n", count, status, (*lscm).hsolver_data.n);
  return status;
}

char coord_pars_corresp(coord_pars *a, coord_pars *b)
{
  // This assumes that (*a).ts == (*b).ts
  for (int i = 0; i < (*a).map.len; i++)
    {
      if ((*a).map.e[i] == (*b).map.e[i]) {}
      else return 0;
    }
  return 1;
}

int lsc_h_minimizer_solve(lsc_h_minimizer *lscm, double tol)
{
  // Store last configuration (excluding topology), measure distance between consecutive minima
  array_int map;
  aarray_int fibers;
  array_int cmobile_map;
  array_int cmobile;
  lsc_heuristic_pars *hpars = &((*lscm).hpars);
  coord_pars *cpars = &((*lscm).cpars);
  coord_pars cpars0;
  transcribe_coord_pars(cpars, &cpars0);
  gsl_vector *x0 = gsl_vector_alloc((*lscm).hsolver_data.n);
  gsl_vector_memcpy(x0, (*lscm).c_data);
  int count = 0;
  double tolsq = tol * tol;
  int status = GSL_CONTINUE;
  while (count < MAX_ITER)
    {
      count += 1;
      int status = lsc_h_minimizer_relax(lscm, tol);
      if (status == GSL_SUCCESS) 
	{
	  // Compute distance to last configuration
	  double distsq = contr_pointset_distsq(x0, &cpars0, (*lscm).hsolver->x, cpars);
	  printf("lsc_h_minimizer_solve (%g): err = %g\n", tol, sqrt(distsq));
	  if (distsq < tolsq)
	    {
	      status = GSL_SUCCESS;
	      break;
	    }
	}
      // Transcribe coordinates (and mapping if necessary)
      if (coord_pars_corresp(&cpars0, cpars)) {}
      else
	{
	  free_coord_pars(&cpars0);
	  gsl_vector_free(x0);
	  transcribe_coord_pars(cpars, &cpars0);
	  x0 = gsl_vector_alloc((*lscm).hsolver_data.n);
	}
      gsl_vector_memcpy(x0, (*lscm).hsolver->x);
      // Reset/expand the solver
      lsc_h_minimizer_expand_clusters(lscm);
      // Perform several iterations (try to refine this step)
      for (int i = 0; i < 5; i++) lsc_h_minimizer_iterate(lscm);
    }
  // Free auxiliary structures
  return status;
}

// NOTE: consider replacing this with a 'mobile' string_length function (and the 'Ls' array
//   in lsc_minimizer with a list of targets for the mobile string length.)
double lsc_string_length(void *ctop_, void *cemb_, coord_pars *cpars, const gsl_vector *x)
{
  edge_wtd_graph *ctop = (edge_wtd_graph *) ctop_;
  //fprintf_nbrlist(&((*ctop).top), stdout);
  array_int *cemb = (array_int *) cemb_;
  //fprintf_array_int(cemb, stdout);
  tack_set *ts = (*cpars).ts;
  double L = 0;
  for (int ti = 0; ti < (*ctop).top.v.len; ti++)
    {
      const double *xci;
      int ci = (*cemb).e[ti];
      int mci = (*cpars).cmobile_map.e[ci];
      if (mci > -1) xci = gsl_vector_const_ptr(x, (*cpars).c_map.e[mci]);
      else
	{
	  int rci = (*cpars).fibers.e[ci].e[0];
	  xci = (double *) (*ts).coords.e[rci];
	}
      for (int ni = 0; ni < (*ctop).top.v.e[ti].len; ni++)
	{
	  int tii = (*ctop).top.v.e[ti].e[ni];
	  if (tii > ti) {}
	  else continue;
	  const double *xcii;
	  int cii = (*cemb).e[tii];
	  int mcii = (*cpars).cmobile_map.e[cii];
	  if (mcii > -1) xcii = gsl_vector_const_ptr(x, (*cpars).c_map.e[mcii]);
	  else
	    {
	      int rcii = (*cpars).fibers.e[cii].e[0];
	      xcii = (double *) (*ts).coords.e[rcii];
	    }
	  double delxsq = euclid_distsq(xci, xcii, (*ts).dim);
	  L += (*ctop).edge_wts.e[ti].e[ni] * sqrt(delxsq);
	}
    }
  return L;
}

double lsc_h_f(const gsl_vector *x, void *pars)
{
  // H = 0.5 * k * sum_{tops} (L[top] - L0[top])^2 + H_,
  //  where H_ = -ext_f . x_i (if ext_f != NULL) or H_ = L[top_i] if ext_f == NULL
  lsc_heuristic_pars *hpars = (lsc_heuristic_pars *) pars;
  coord_pars *cpars = (*hpars).cpars;
  double H_ = 0;
  linked_sc *lsc = (*hpars).lsc;
  tack_set *ts = (*lsc).ts;
  if ((*hpars).ext_f == NULL) H_ = lsc_string_length((*(*hpars).ctops).e[(*hpars).i], (*(*hpars).cembs).e[(*hpars).i], cpars, x);
  else
    {
      const double *xi;
      int ci = (*cpars).map.e[(*hpars).i];
      int mci = (*cpars).cmobile_map.e[ci];
      if (mci > -1) xi = gsl_vector_const_ptr(x, (*cpars).c_map.e[mci]);
      else
	{
	  printf("Something weird happened! (this message should not appear!\n");
	  int ri = (*cpars).fibers.e[ci].e[0];
	  xi = (double *) (*ts).coords.e[ri];
	}
      for (int di = 0; di < (*ts).dim; di++) H_ -= (*hpars).ext_f[di] * xi[di];
    }
  // Compute contributions from non-minimized string segments
  if ((*hpars).Ls != NULL)
    {
      double H_1 = 0;
      for (int ti = 0; ti < ((*hpars).ctops)->len; ti++)
	{
	  if ((*hpars).ext_f != NULL || (*hpars).i != ti) {}
	  else continue;
	  double L_ti = lsc_string_length((*(*hpars).ctops).e[ti], (*(*hpars).cembs).e[ti], cpars, x);
	  L_ti -= (*hpars).Ls[ti];
	  H_1 += L_ti * L_ti;
	}
      H_1 *= 0.5 * (*hpars).k;
      H_ += H_1;
    }
  return H_;
}

void lsc_h_df_aux(const gsl_vector *x, lsc_heuristic_pars *hpars, coord_pars *cpars, tack_set *ts, double *prefactors, gsl_vector *df)
{
  for (int mcii = 0; mcii < (*cpars).cmobile.len; mcii++)
    {
      int cii = (*cpars).cmobile.e[mcii];
      const double *xcii = gsl_vector_const_ptr(x, (*cpars).c_map.e[mcii]);
      double *fcii = gsl_vector_ptr(df, (*cpars).c_map.e[mcii]);
      for (int ti = 0; ti < (*hpars).ctops->len; ti++)
	{
	  int tii = (*(*hpars).cembs_map).e[cii].e[ti];
	  if (tii > -1) {}
	  else continue;
	  edge_wtd_graph *ctop = (edge_wtd_graph *) (*hpars).ctops->e[ti];
	  array_int *cemb = (array_int *) (*hpars).cembs->e[ti];
	  double prefactor = ((*hpars).ext_f != NULL || ti != (*hpars).i ? prefactors[ti] : 1);
	  for (int ntii = 0; ntii < (*ctop).top.v.e[tii].len; ntii++)
	    {
	      int tiii = (*ctop).top.v.e[tii].e[ntii];
	      int ciii = (*cemb).e[tiii];
	      int mciii = (*cpars).cmobile_map.e[ciii];
	      if (mciii > mcii) continue;
	      const double *xciii;
	      double *fciii = NULL;
	      if (mciii > -1)
		{
		  xciii = gsl_vector_const_ptr(x, (*cpars).c_map.e[mciii]);
		  fciii = gsl_vector_ptr(df, (*cpars).c_map.e[mciii]);
		}
	      else
		{
		  int rciii = (*cpars).fibers.e[ciii].e[0];
		  xciii = (double *) (*ts).coords.e[rciii];
		}
	      double disp[(*ts).dim];
	      double dispsq = 0;
	      for (int di = 0; di < (*ts).dim; di++)
		{
		  disp[di] = xcii[di] - xciii[di];
		  dispsq += disp[di] * disp[di];
		}
	      if (dispsq > 0) {}
	      else continue;
	      dispsq = prefactor * (*ctop).edge_wts.e[tii].e[ntii] / sqrt(dispsq);
	      for (int di = 0; di < (*ts).dim; di++)
		{
		  disp[di] *= dispsq;
		  fcii[di] += disp[di];
		}
	      if (fciii != NULL)
		{
		  for (int di = 0; di < (*ts).dim; di++) fciii[di] -= disp[di];
		}
	    }
	}
    }
    // Add gradient of H_
    if ((*hpars).ext_f != NULL)
    {
      int ci = (*cpars).map.e[(*hpars).i];
      int mci = (*cpars).cmobile_map.e[ci];
      if (mci > -1) {}
      else
	{
	  printf("Something weird happened! (This message should not appear!)\n");
	  exit(EXIT_FAILURE);
	}
      double *df_ci = gsl_vector_ptr(df, (*cpars).c_map.e[mci]);
      for (int di = 0; di < (*ts).dim; di++) df_ci[di] -= (*hpars).ext_f[di];
    }
}

void lsc_h_df(const gsl_vector *x, void *pars, gsl_vector *df)
{
  lsc_heuristic_pars *hpars = (lsc_heuristic_pars *) pars;
  coord_pars *cpars = (*hpars).cpars;
  linked_sc *lsc = (*hpars).lsc;
  tack_set *ts = (*lsc).ts;
  // Compute gradient of (approximate) length constraint
  //  dH = k * sum_{tops} (L[top] - L0) * dL[top] 
  gsl_vector_set_zero(df);
  double prefactors[(*hpars).ctops->len];
  if ((*hpars).Ls != NULL)
    {
      for (int ti = 0; ti < (*hpars).ctops->len; ti++)
	{
	  if ((*hpars).ext_f != NULL || ti != (*hpars).i) prefactors[ti] = (*hpars).k * (lsc_string_length((*(*hpars).ctops).e[ti], (*(*hpars).cembs).e[ti], cpars, x) - (*hpars).Ls[ti]);
	}
    }
  lsc_h_df_aux(x, hpars, cpars, ts, &(prefactors[0]), df);
}

void lsc_h_fdf(const gsl_vector *x, void *pars, double *f, gsl_vector *df)
{
  lsc_heuristic_pars *hpars = (lsc_heuristic_pars *) pars;
  linked_sc *lsc = (*hpars).lsc;
  tack_set *ts = (*lsc).ts;
  coord_pars *cpars = (*hpars).cpars;
  (*f) = 0;
  gsl_vector_set_zero(df);
  double prefactors[(*hpars).ctops->len];
  if ((*hpars).Ls != NULL)
    {
      for (int ti = 0; ti < (*hpars).ctops->len; ti++)
	{
	  //(lsc_string_length((*(*hpars).ctops).e[ti], (*(*hpars).cembs).e[ti], cpars, x) - (*hpars).Ls[ti]);
	  double incr =  lsc_string_length((*(*hpars).ctops).e[ti], (*(*hpars).cembs).e[ti], cpars, x) - (*hpars).Ls[ti];
	  if ((*hpars).ext_f != NULL || (*hpars).i != ti)
	    {
	      prefactors[ti] = (*hpars).k * incr;
	      (*f) += prefactors[ti] * incr;
	    }
	}
      (*f) *= 0.5;
    }
  if ((*hpars).ext_f != NULL)
    {
      int ci = (*cpars).map.e[(*hpars).i];
      int mci = (*cpars).cmobile_map.e[ci];
      const double *xi = gsl_vector_const_ptr(x, (*cpars).c_map.e[mci]);
      for (int di = 0; di < (*ts).dim; di++)
	{
	  (*f) -= (*hpars).ext_f[di] * xi[di];
	}
    }
  else (*f) += lsc_string_length((*(*hpars).ctops).e[(*hpars).i], (*(*hpars).cembs).e[(*hpars).i], cpars, x);
  lsc_h_df_aux(x, hpars, cpars, ts, &(prefactors[0]), df);
}

void lsc_L_solver_init(lsc_L_solver *lscs, lsc_h_minimizer *lscm)
{
  // TBD
}

void free_lsc_L_solver(lsc_L_solver *lscs)
{
  // TBD
}

void lsc_L_solver_iterate(lsc_L_solver *lscs)
{
  //...
}

void lsc_L_solver_relax(lsc_L_solver *lscs, double tol)
{

}
void lsc_L_solver_solve(lsc_L_solver *lscs, double tol)
{

}

int lsc_L_f(const gsl_vector *x, void *pars, gsl_vector *f)
{
}
int lsc_L_df(const gsl_vector *x, void *pars, gsl_matrix *J)
{

}
int lsc_L_fdf(const gsl_vector *x, void *pars, gsl_vector *f, gsl_matrix *J)
{

}
