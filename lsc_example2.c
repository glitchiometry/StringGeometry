#include "sconfig.h"
#include "gen_configs.h"

double tol = 1e-7;

int N_samples = 40;

int main()
{
  array_voidstar coords;
  /*
     Generate "tack set" in the pattern:
     .                       1
       o  o                    3  0
     .        with indexing  2
   (dots and 'o''s represent fixed and mobile vertices respectively)
  */
  init_coords_Y(&coords);
  tack_set ts;
  tack_set_init(&ts, 2);
  add_mbl_pt_tack_set(&ts, (double *) coords.e[0]);
  for (int i = 1; i < 3; i++) add_fxd_pt_tack_set(&ts, (double *) coords.e[i]);
  add_mbl_pt_tack_set(&ts, (double *) coords.e[3]);
  // Define the asssociated linked string configuration
  linked_sc lsc;
  linked_sc_init(&lsc, &ts);
  // Add topologies and embeddings to the linked string configuration 
  // Some care is needed here to ensure
  //    (1) All vertices in the 'tack set' are covered and connected to each other
  //    (2) The topologies connect at at least one mobile vertex.
  //    (3) When an external force is applied, either (i) topologies connect at two or more mobile vertices,
  //       or (ii) the force is applied at a mobile vertex on only one of the topologies (which connects to
  //       the others on at least one additional mobile site/vertex.)
  //    (4) That embeddings and topologies/graphs are defined consistently.
  /* 
     Define 'V' graphs with the following embeddings:
    1: .      2:  
       \        '  
        o--o     o--o
     .          /
               '
   */
  edge_wtd_graph top1;
  array_int emb1;
  init_top_elbow(&top1);
  array_int_init(&emb1, 3);
  emb1.len = 3;
  emb1.e[0] = 0;
  emb1.e[1] = 3;
  emb1.e[2] = 1;

  add2linked_sc(&lsc, &top1, &emb1);
  array_int emb2;
  edge_wtd_graph top2;
  init_top_elbow(&top2);
  array_int_init(&emb2, 3);
  emb2.e[0] = 0;
  emb2.e[1] = 3;
  emb2.e[2] = 2;
  emb2.len = 3;
  add2linked_sc(&lsc, &top2, &emb2);
  // The external force will be applied at vertex '1' (the top/center site.)
  int ext_f_site = 0;
  double Ls[2] = {1.2, 1.2};
  double dtheta = (2 * M_PI) / N_samples;
  double theta = 0;
  FILE *ofile = fopen("lsc_example2.dat", "w");
  for (int si = 0; si < N_samples; si++)
    {
      lsc_h_minimizer lscm;
      double ext_f[2] = {cos(theta), sin(theta)};
      lsc_h_minimizer_init(&lscm, &lsc, ext_f_site, &ext_f[0], &Ls[0], 1.0);
      fprintf(ofile, "%g ", theta);
      int status = lsc_h_minimizer_solve(&lscm, tol);
      for (int i = 0; i < ts.coords.len; i++)
	{
	  int ci = lscm.cpars.map.e[i];
	  int mci = lscm.cpars.cmobile_map.e[ci];
	  double *xci;
	  if (mci > -1) xci = gsl_vector_ptr(lscm.hsolver->x, lscm.cpars.c_map.e[mci]);
	  else
	    {
	      int ri = lscm.cpars.fibers.e[ci].e[0];
	      xci = (double *) ts.coords.e[ri];
	    }
	  for (int di = 0; di < ts.dim; di++) fprintf(ofile, "%g ", xci[di]);
	}
      theta += dtheta;
      fprintf(ofile, "\n");
      free_lsc_h_minimizer(&lscm);
    }
  fclose(ofile);
  free_linked_sc_shallow(&lsc);
  free_array_int(&emb1);
  free_array_int(&emb2);
  free_edge_wtd_graph(&top1);
  free_edge_wtd_graph(&top2);

  free_tack_set(&ts, NULL);
  free_array_voidstar(&coords, free_elem_triv);
  return 0;
}
