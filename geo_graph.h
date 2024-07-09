#ifndef GEO_GRAPH_H
#define GEO_GRAPH_H
#include "basics.h"

typedef struct
{
	nbrlist top;
	aarray_double pos;
	int dim;
	aarray_int wts;
} geo_graph;

void geo_graph_init(geo_graph *gg, int N, int dim);
void free_geo_graph(geo_graph *gg);
void add_point_geo_graph(geo_graph *gg, double *x);
void remove_point_geo_graph(geo_graph *gg, int i);
void add_edge_geo_graph(geo_graph *gg, int i, int j, int wt);
int add_edge_geo_graph_safe(geo_graph *gg, int i, int j, int wt);
void add_edge_wt_geo_graph(geo_graph *gg, int i, int j, int wt);
void remove_edge_geo_graph(geo_graph *gg, int i, int ni);
char geo_graph_check_consistent_wts(geo_graph *gg);
char geo_graph_check_consistent_wts_edge(geo_graph *gg, int i, int ni);


double geo_graph_dist(geo_graph *gg, int i, int j);
double geo_graph_distsq(geo_graph *gg, int i, int j);
void geo_graph_disp(geo_graph *gg, int i, int j, double *xjmxi);
void fprintf_geo_graph(geo_graph *gg, FILE *ofile);
#endif
