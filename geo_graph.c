#include "geo_graph.h"

void geo_graph_init(geo_graph *gg, int N, int dim)
{
	nbrlist_init_precise(&((*gg).top), N);
	aarray_double_init_precise(&((*gg).pos), N, dim);
	aarray_int_init(&((*gg).wts), N);
	(*gg).dim = dim;
}

void free_geo_graph(geo_graph *gg)
{
	free_aarray_double(&((*gg).pos));
	free_nbrlist(&((*gg).top));
	free_aarray_int(&((*gg).wts));
}

void add_point_geo_graph(geo_graph *gg, double *x)
{
	array_double ax;
	array_double_init(&ax, (*gg).dim);
	for (int i = 0; i < (*gg).dim; i++) ax.e[i] = x[i];
	add2aarray_double(&((*gg).pos), ax);
	extend_nbrlist(&((*gg).top));
	extend_aarray_int(&((*gg).wts));
}

void remove_point_geo_graph(geo_graph *gg, int i)
{
	remove_aarray_double(&((*gg).pos), i);
	remove_vertex_nbrlist(&((*gg).top), i);
	remove_aarray_int(&((*gg).wts), i);
}

int add_edge_geo_graph_safe(geo_graph *gg, int i, int j, int wt)
{
	int ni = add_edge_nbrlist_safe(&((*gg).top), i, j);
	if (ni > -1) 
	{
		return ni;
	}
	else
	{
		add2array_int(&((*gg).wts.e[i]), wt);
		add2array_int(&((*gg).wts.e[j]), wt);
		return -1;
	}
}

void add_edge_geo_graph(geo_graph *gg, int i, int j, int wt)
{
	add_edge_nbrlist(&((*gg).top), i, j);
	add2array_int(&((*gg).wts.e[i]), wt);
	add2array_int(&((*gg).wts.e[j]), wt);
}

void add_edge_wt_geo_graph(geo_graph *gg, int i, int j, int wt)
{
  int ni = add_edge_geo_graph_safe(gg, i, j, wt);
  if (ni > -1)
    {
      int i_i_j = (*gg).top.i_of.e[i].e[ni];
      (*gg).wts.e[i].e[ni] += wt;
      (*gg).wts.e[j].e[i_i_j] += wt;
    }
}

void remove_edge_geo_graph(geo_graph *gg, int i, int ni)
{
	int j = (*gg).top.v.e[i].e[ni];
	int i_ij = (*gg).top.i_of.e[i].e[ni];
	remove_array_int(&((*gg).wts.e[i]), ni);
	remove_array_int(&((*gg).wts.e[j]), i_ij);
	remove_edge_nbrlist(&((*gg).top), i, ni);
}

double geo_graph_dist(geo_graph *gg, int i, int j)
{
	return sqrt(geo_graph_distsq(gg, i, j));
}

double geo_graph_distsq(geo_graph *gg, int i, int j)
{
	double delsq = 0;
	for (int di = 0; di < (*gg).dim; di++)
	{
		double delx = (*gg).pos.e[i].e[di] - (*gg).pos.e[j].e[di];
		delsq += delx * delx;
	}
	return delsq;
}

void geo_graph_disp(geo_graph *gg, int i, int j, double *xjmxi)
{
	for (int di = 0; di < (*gg).dim; di++)
	{
		xjmxi[di] = (*gg).pos.e[j].e[di] - (*gg).pos.e[i].e[di];
	}
}

double geo_graph_disp_distsq(geo_graph *gg, int i, int j, double *xjmxi)
{
	double distsq = 0;
	for (int di = 0; di < (*gg).dim; di++)
	{
		xjmxi[di] = (*gg).pos.e[j].e[di] - (*gg).pos.e[i].e[di];
		distsq += xjmxi[di] * xjmxi[di];
	}
	return distsq;
}

char geo_graph_check_consistent_wts_edge(geo_graph *gg, int i, int ni)
{
	int j = (*gg).top.v.e[i].e[ni];
	int i_ij = (*gg).top.i_of.e[i].e[ni];
	if ((*gg).wts.e[i].e[ni] == (*gg).wts.e[j].e[i_ij]) return 1;
	else return 0;
}

char geo_graph_check_consistent_wts(geo_graph *gg)
{
	for (int i = 0; i < (*gg).top.v.len; i++)
	{
		for (int ni = 0; ni < (*gg).top.v.e[i].len; ni++) 
		{
			int j = (*gg).top.v.e[i].e[ni];
			if (j < i) continue;
			geo_graph_check_consistent_wts_edge(gg, i, ni);	
		}
	}
	return 1;
}

void fprintf_geo_graph(geo_graph *gg, FILE *ofile)
{
	if (ofile != NULL)
	{
		printf("%d %d\n", (*gg).top.v.len, (*gg).dim);
		for (int i = 0; i < (*gg).top.v.len; i++)
		{
			for (int di = 0; di < (*gg).dim; di++)
			{
				printf("%g ", (*gg).pos.e[i].e[di]);
			}
			for (int ni = 0; ni < (*gg).top.v.e[i].len; ni++)
			{
				printf("%d ", (*gg).top.v.e[i].e[ni]);
			}
			printf("\n");
		}
	}
}
