#include "gen_configs.h"

void string_config_init_H(string_config *sc, array_voidstar *coords, double width, double height)
{
  string_config_init(sc, 6, 2);
  array_voidstar_init(coords, 6);
  double *x[6];
  for (int i = 0; i < 6; i++)
    {
      double *xi = (double *) malloc(sizeof(double) << 1);
      add2array_voidstar(coords, xi);
      set_pos_string_config(sc, i, xi);
      x[i] = xi;
    }
  x[2][0] = x[3][0] = x[0][0] = 0;
  x[2][1] = x[4][1] = height;
  x[3][1] = x[5][1] = 0;
  x[5][0] = x[4][0] = x[1][0] = width;
  x[0][1] = x[1][1] = 0.5 * height;
  add_edge_string_config(sc, 0, 2, 2);
  add_edge_string_config(sc, 0, 3, 2);
  add_edge_string_config(sc, 0, 1, 2);
  add_edge_string_config(sc, 1, 4, 2);
  add_edge_string_config(sc, 1, 5, 2);
  for (int i = 2; i < 6; i++) fix_point_string_config(sc, i);
}

void string_config_init_Y(string_config *sc, array_voidstar *coords)
{
  string_config_init(sc, 4, 2);
  array_voidstar_init(coords, 4);
  for (int i = 0; i < 4; i++)
    {
      double *xi = (double *) calloc(2, sizeof(double));
      add2array_voidstar(coords, xi);
      set_pos_string_config(sc, i, xi);
    }
  double dtheta = (2 * M_PI) / 3;
  double theta = 0;
  for (int i = 0; i < 3; i++)
    {
      add_edge_string_config(sc, i, 3, 2);
      fix_point_string_config(sc, i);
      double *xi = (double *) (*coords).e[i];
      xi[0] = cos(theta);
      xi[1] = sin(theta);
      theta += dtheta;
    }
  double *xm = (double *) (*coords).e[3];
  xm[0] = 0.3;
  xm[1] = 0.3;
}

void string_config_init_10(string_config *sc, array_voidstar *coords)
{
  string_config_init(sc, 2, 2);
  array_voidstar_init(coords, 2);
  for (int i = 0; i < 2; i++)
    {
      double *xi = (double *) calloc(2, sizeof(double));
      add2array_voidstar(coords, xi);
      set_pos_string_config(sc, i, xi);
      xi[0] = (double) i;
      xi[1] = 0;
    }
  fix_point_string_config(sc, 0);
  add_edge_string_config(sc, 0, 1, 2);
}

void string_config_init_X(string_config *sc, array_voidstar *coords, double width, double height)
{
  string_config_init(sc, 5, 2);
  array_voidstar_init(coords, 5);
  double *x[5];
  for (int i = 0; i < 5; i++)
    {
      x[i] = (double *) malloc(sizeof(double)<<1);
      add2array_voidstar(coords, x[i]);
      set_pos_string_config(sc, i, x[i]);
    }
  for (int i = 0; i < 4; i++)
    {
      fix_point_string_config(sc, i);
      add_edge_string_config(sc, i, 4, 2);
    }
  x[0][0] = x[3][0] = width;
  x[0][1] = x[1][1] = height;
  x[1][0] = x[2][0] = x[2][1] = x[3][1] = 0;
}

void string_config_init_101(string_config *sc, array_voidstar *coords)
{
  string_config_init(sc, 3, 2);
  array_voidstar_init(coords, 3);
  for (int i = 0; i < 3; i++)
    {
      double *xi = (double *) calloc(2, sizeof(double));
      add2array_voidstar(coords, xi);
      set_pos_string_config(sc, i, xi);
      xi[0] = i * 0.5;
      xi[1] = (i & 1) * 0.8;
    }
  fix_point_string_config(sc, 0);
  fix_point_string_config(sc, 2);
  add_edge_string_config(sc, 0, 1, 1);
  add_edge_string_config(sc, 2, 1, 1);
}

void string_config_init_random_loop(string_config *sc, array_voidstar *coords, int N_pts, int dim, double side_len)
{
  array_voidstar_init(coords, N_pts);
  string_config_init(sc, N_pts, dim);
  for (int i = 0; i < N_pts; i++)
    {
      double *xi = (double *) calloc(dim, sizeof(double));
      set_pos_string_config(sc, i, xi);
      add2array_voidstar(coords, xi);
      for (int di = 0; di < dim; di++) xi[di] = rnd() * side_len;
    }
  union_find uf;
  union_find_init(&uf, N_pts);
  int i0 = rand() % N_pts;
  int i = i0;
  while (1)
    {
      int j = rand() % N_pts;
      if (j != i) {}
      else continue;
      int ni = add_edge_string_config(sc, i, j, 1);
      if (ni == -1)
	{
	  union_find_union(&uf, i, j);
	  // Check if component of 'i' has size N_pts
	  int comp_size = union_find_order(&uf, i);
	  if (comp_size < N_pts) {}
	  else break;
	}
      i = j;
    }
  if (i != i0) // This case should always hold
    {
      add_edge_string_config(sc, i, i0, 1);
    }
  free_union_find(&uf);
}

void string_config_init_random_loop_stable(string_config *sc, array_voidstar *coords, int N, int dim, double side_len)
{
  string_config_init_random_loop(sc, coords, N, dim, side_len);
  if (N > 3) {}
  else return;
  array_int unstable;
  array_int_init(&unstable, 1);
  while (1) // This should finish after one iteration
    {
      // Find vertices with less than three neighbors
      unstable.len = 0;
      for (int i = 0; i < N; i++)
	{
	  if ((*sc).top.v.e[i].len < 3)
	    {
	      add2array_int(&unstable, i);
	    }
	}
      if (unstable.len > 0) {}
      else break;
      for (int i = 0; i < unstable.len; i++)
	{
	  int ii = unstable.e[i];
	  int jj;
	  do
	    {
	      jj = rand() % N;
	    }
	  while (jj == ii);
	  add_edge_string_config(sc, ii, jj, 2);
	}
    }
  free_array_int(&unstable);
}

void string_config_init_random(string_config *sc, array_voidstar *coords, int N, int dim, int max_wt)
{
  string_config_init(sc, N, dim);
  array_voidstar_init(coords, N);
  for (int i = 0; i < N; i++)
    {
      double *xi = (double *) calloc(dim, sizeof(double));
      for (int di = 0; di < dim; di++) xi[di] = rnd();
      set_pos_string_config(sc, i, xi);
      add2array_voidstar(coords, xi);
    }
  union_find uf;
  union_find_init(&uf, N);
  while (1)
    {
      int i = rand() % N;
      int j = rand() % N;
      if (j == i) continue;
      int wt = 1 + (rand() % max_wt);
      int ni = add_edge_string_config(sc, i, j, wt);
      union_find_union(&uf, i, j);
      int ci;
      int dci = union_find_find(&uf, i, &ci);
      if (uf.cluster_size.e[ci] == N) break;
    }
  free_union_find(&uf);
}

void init_coords_H(array_voidstar *coords, double width, double height)
{
  array_voidstar_init(coords, 6);
  double *x[6];
  for (int i = 0; i < 6; i++)
    {
      x[i] = (double *) malloc(sizeof(double) << 1);
      add2array_voidstar(coords, x[i]);
    }
  x[2][0] = x[3][0] = x[0][0] = 0;
  x[2][1] = x[4][1] = height;
  x[3][1] = x[5][1] = 0;
  x[5][0] = x[4][0] = x[1][0] = width;
  x[0][1] = x[1][1] = 0.5 * height;
}

void init_coords_Y(array_voidstar *coords)
{
  array_voidstar_init(coords, 4);
  for (int i = 0; i < 4; i++)
    {
      double *xi = (double *) calloc(2, sizeof(double));
      add2array_voidstar(coords, xi);
    }
  double dtheta = (2 * M_PI) / 3;
  double theta = 0;
  for (int i = 0; i < 3; i++)
    {
      double *xi = (double *) (*coords).e[i];
      xi[0] = cos(theta);
      xi[1] = sin(theta);
      theta += dtheta;
    }
  double *xm = (double *) (*coords).e[3];
  xm[0] = 0.3;
  xm[1] = 0.3;
}

void init_coords_X(array_voidstar *coords, double width, double height)
{
  array_voidstar_init(coords, 5);
  double *x[5];
  for (int i = 0; i < 5; i++)
    {
      x[i] = (double *) malloc(sizeof(double)<<1);
      add2array_voidstar(coords, x[i]);
    }
  x[0][0] = x[3][0] = width;
  x[0][1] = x[1][1] = height;
  x[1][0] = x[2][0] = x[2][1] = x[3][1] = 0;
}

void init_coords_V(array_voidstar *coords)
{
  array_voidstar_init(coords, 3);
  double *x[3];
  for (int i = 0; i < 3; i++)
    {
      x[i] = (double *) malloc(sizeof(double)<<1);
      add2array_voidstar(coords, x[i]);
    }
  x[0][0] = x[0][1] = x[2][1] = 0;
  x[2][0] = 1;
  x[1][1] = -(x[1][0] = 0.5);
}

void init_coords_grid(array_voidstar *coords, int N, int M, double a)
{
  int N_total = N * M;
  double *x[N][M];
  array_voidstar_init(coords, N_total);
  (*coords).len = N_total;
  int count = 0;
  for (int i = 0; i < N; i++)
    {
      for (int ii = 0; ii < M; ii++)
	{
	  x[i][ii] = malloc(sizeof(double) << 1);
	  x[i][ii][0] = i * a;
	  x[i][ii][1] = ii * a;
	  (*coords).e[count] = x[i][ii];
	  count += 1;
	}
    }
}

void init_coords_hex(array_voidstar *coords, double R, double a)
{
  double Rsq = R * R;
  double v[3][2] = {{a, 0}, {-0.5 * a, sqrt(0.75) * a}, {-1, -1}};
  v[2][0] = v[1][0];
  v[2][1] = -v[1][1];
  int v_[3][2] = {{1, 0}, {-1, 1}, {0, -1}};
  array_int bdry;
  array_int_init(&bdry, 1);
  hash_table_int_str ht; // 3 * 0.25 * a^2 * sqrt(3)
  double aux_radius = R / a;
  int int_radius = (int) (aux_radius) + 1;
  int size_est = (int) ((M_PI * aux_radius * aux_radius) / sqrt(1.6875));
  array_voidstar_init(coords, size_est);
  hash_table_int_str_init(&ht, size_est + 1, int_radius, BASICS_H_SELF);
  add2array_int(&bdry, 0);
  double *x0 = malloc(sizeof(double) << 1);
  x0[0] = x0[1] = 0;
  add2array_voidstar(coords, x0);
  array_voidstar coords_;
  array_char parity;
  array_voidstar_init(&coords_, 1);
  array_char_init(&parity, 1);
  int *addr = malloc(sizeof(int)<<1);
  addr[0] = addr[1] = 0;
  add2array_voidstar(&coords_, addr);
  add2array_char(&parity, 0);
  add2hash_table_int_str(&ht, addr, 2, NULL);
  while (bdry.len > 0)
    {
      bdry.len -= 1;
      int i = bdry.e[bdry.len];
      int *addr_i = (int *) coords_.e[i];
      int aux_addr[2];
      int phase = parity.e[i] == 0 ? 1 : -1;
      for (int ii = 0; ii < 3; ii++)
	{
	  int aux_addr[2];
	  aux_addr[0] = addr_i[0] + phase * v_[ii][0];
	  aux_addr[1] = addr_i[1] + phase * v_[ii][1];
	  int elem_addr;
	  char present = query_hash_table_int_str(&ht, &aux_addr[0], 2, &elem_addr, NULL);
	  if (present) {}
	  else
	    {
	      // Check if the proposed site is within radius 'R' of the origin
	      double *x_i = (double *) coords->e[i];
	      double x_aux[2] = {x_i[0] + phase * v[ii][0], x_i[1] + phase * v[ii][1]};
	      if ((x_aux[0] * x_aux[0] + x_aux[1] * x_aux[1]) < Rsq)
		{
		  add2array_int(&bdry, coords_.len);
		  add2array_char(&parity, !parity.e[i]);
		  int *addr_aux = malloc(sizeof(int)<<1);
		  addr_aux[0] = aux_addr[0];
		  addr_aux[1] = aux_addr[1];
		  add2array_voidstar(&coords_, addr_aux);
		  double *x_ni = malloc(sizeof(double)<<1);
		  x_ni[0] = x_aux[0];
		  x_ni[1] = x_aux[1];
		  add2array_voidstar(coords, x_ni);
		  add2hash_table_int_str(&ht, addr_aux, 2, NULL);
		}
	    }
	}
    }
  free_array_voidstar(&coords_, free_elem_triv);
  free_array_char(&parity);
  free_array_int(&bdry);
  free_hash_table_int_str(&ht, NULL);
}

void init_coords_loop(array_voidstar *coords, int N)
{
  double dtheta = (2 * M_PI ) / N;
  double theta = 0;
  array_voidstar_init(coords, N);
  for (int i = 0; i < N; i++)
    {
      double *xi = malloc(sizeof(double)<<1);
      xi[0] = cos(theta);
      xi[1] = sin(theta);
      add2array_voidstar(coords, xi);
      theta += dtheta;
    }
}

void init_coords_A(array_voidstar *coords, double width, double height)
{
  array_voidstar_init(coords, 7);
  double *x[7];
  for (int i = 0; i < 7; i++)
    {
      x[i] = (double *) malloc(sizeof(double) << 1);
      add2array_voidstar(coords, x[i]);
    }
  x[3][0] = x[4][0] = x[0][0] = 0;
  x[3][1] = x[5][1] = x[1][1] = height;
  x[4][1] = x[6][1] = 0;
  x[6][0] = x[5][0] = x[2][0] = width;
  x[0][1] = x[2][1] = 0.5 * height;
  x[1][0] = 0.5 * width;
}

// Topologies
void init_top_H(edge_wtd_graph *top)
{
  edge_wtd_graph_init(top, 6);
  for (int i = 0; i < 6; i++) extend_edge_wtd_graph(top);
  for (int i = 2; i < 4; i++) add_edge_edge_wtd_graph(top, 0, i, 2);
  for (int i = 4; i < 6; i++) add_edge_edge_wtd_graph(top, 1, i, 2);
  add_edge_edge_wtd_graph(top, 0, 1, 2);
}

void init_top_elbow(edge_wtd_graph *top)
{
  edge_wtd_graph_init(top, 3);
  for (int i = 0; i < 3; i++) extend_edge_wtd_graph(top);
  int im1 = 0;
  for (int i = 1; i < 3; i++)
    {
      add_edge_edge_wtd_graph(top, im1, i, 1);
      im1 = i;
    }
}

void init_top_Y(edge_wtd_graph *top)
{
  edge_wtd_graph_init(top, 4);
  for (int i = 0; i < 4; i++) extend_edge_wtd_graph(top);
  for (int i = 0; i < 3; i++) add_edge_edge_wtd_graph(top, 3, i, 2);
}

void init_top_X(edge_wtd_graph *top)
{
  edge_wtd_graph_init(top, 5);
  for (int i = 0; i < 5; i++) extend_edge_wtd_graph(top);
  for (int i = 0; i < 4; i++) add_edge_edge_wtd_graph(top, i, 4, 2);
}

void init_top_dimer(edge_wtd_graph *top)
{
  edge_wtd_graph_init(top, 2);
  extend_edge_wtd_graph(top);
  extend_edge_wtd_graph(top);
  add_edge_edge_wtd_graph(top, 0, 1, 2);
}

void init_top_loop(edge_wtd_graph *top, int N)
{
  edge_wtd_graph_init(top, N);
  for (int i = 0; i < N; i++) extend_edge_wtd_graph(top);
  int im1 = 0;
  for (int i = 1; i < N; i++)
    {
      add_edge_edge_wtd_graph(top, i, im1, 1);
    }
  add_edge_edge_wtd_graph(top, 0, im1, 1);
}
