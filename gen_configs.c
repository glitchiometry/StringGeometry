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
