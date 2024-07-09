#include "random_geo_graph.h"

double rnd()
{
	return ((double) rand()) / RAND_MAX;
}

void random_connected_geo_graph(geo_graph *g, int N_pts, int dim, double side_len)
{
	geo_graph_init(g, N_pts, dim);
	for (int i = 0; i < N_pts; i++)
	{
		double x_i[dim];
		for (int di = 0; di < dim; di++) x_i[di] = rnd() * side_len;
		add_point_geo_graph(g, &(x_i[0]));
	}
	if (N_pts > 1) {}
	else return;
	union_find uf;
	union_find_init(&uf, N_pts);
	while (1)
	{
		int i = rand() % N_pts;
		int j = rand() % N_pts;
		if (j != i) {}
		else continue;
		int ni = add_edge_geo_graph_safe(g, i, j, 1);
		union_find_union(&uf, i, j);
		// Check if component of 'i' has size N_pts
		int comp_size = union_find_order(&uf, i);
		if (comp_size < N_pts) {}
		else break;
	}
	free_union_find(&uf);
}

void random_loop_geo_graph(geo_graph *g, int N_pts, int dim, double side_len)
{
  	geo_graph_init(g, N_pts, dim);
	for (int i = 0; i < N_pts; i++)
	{
		double x_i[dim];
		for (int di = 0; di < dim; di++) x_i[di] = rnd() * side_len;
		add_point_geo_graph(g, &(x_i[0]));
	}
	if (N_pts > 1) {}
	else return;
	union_find uf;
	union_find_init(&uf, N_pts);
	int i = rand() % N_pts;
	while (1)
	{
		int j = rand() % N_pts;
		if (j != i) {}
		else continue;
		int ni = add_edge_wt_geo_graph(g, i, j, 1);
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
	free_union_find(&uf);
}
