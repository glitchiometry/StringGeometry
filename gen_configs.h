#include "sconfig.h"

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

