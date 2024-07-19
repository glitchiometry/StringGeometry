#include "sconfig.h"

// Coordinates
void init_coords_H(array_voidstar *coords, double width, double height);
void init_coords_Y(array_voidstar *coords);
void init_coords_X(array_voidstar *coords, double width, double height);
void init_coords_V(array_voidstar *coords);
void init_coords_grid(array_voidstar *coords, int N, int M, double a);
void init_coords_hex(array_voidstar *coords, double R, double a);
void init_coords_loop(array_voidstar *coords, int N);
void init_coords_A(array_voidstar *coords, double width, double height);

// Topologies
void init_top_H(edge_wtd_graph *top);
void init_top_elbow(edge_wtd_graph *top);
void init_top_Y(edge_wtd_graph *top);
void init_top_X(edge_wtd_graph *top);
void init_top_dimer(edge_wtd_graph *top);
void init_top_loop(edge_wtd_graph *top, int N);

// Input: uninitialized string config 
void string_config_init_10(string_config *sc, array_voidstar *coords);
void string_config_init_Y(string_config *sc, array_voidstar *coords);
void string_config_init_101(string_config *sc, array_voidstar *coords);
void string_config_init_H(string_config *sc, array_voidstar *coords, double width, double height);
void string_config_init_X(string_config *sc, array_voidstar *coords, double width, double height);
void string_config_init_random_loop(string_config *sc, array_voidstar *coords, int N_pts, int dim, double side_len);
void string_config_init_random_loop_stable(string_config *sc, array_voidstar *coords, int N, int dim, double side_len);
void string_config_init_random(string_config *sc, array_voidstar *coords, int N, int dim, int max_wt);
void fprintf_string_config(string_config *sc, char *dirname);

