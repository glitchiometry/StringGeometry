String Geometry: a repository for "string and thumbtack" constructions (in dimension > 1)

Requirements: 
  - GNU Scientific Library (GSL: https://www.gnu.org/software/gsl/)
    
  from whw_clib:
  - basics.c/h
  - system.c/h (only for saving configurations)

TO DO:
- Change sc_minimizer_mm (multimin) and sc_Lagrange_solver to track indirectly frozen vertices (i.e. vertices whose neighbors are all frozen)
  - This might improve numerical stability of gsl_multiroot_fdf_solvers.
- Test performance and precision of merging/unmerging criteria in sc_minimizer_mm.
- Devise multi-scale variants of current algorithms to resolve near-singular vertex arrangements (and facilitate potential applications in mesoscopic design.)
- Devise a system for handling "second order" singular configurations (i.e. where the Jacobian is singular at a root.)
- Implement an 'adiabatic' variant of the constant length (Lagrange multiplier) solver, to facilitate computing samples for representing (approximate) constant length hypersurfaces
- Implement an interpolant to represent (in conjunction with the Lagrange solver) constant length surfaces
- Implement a generic Lagrange multiplier solver for coupled disjoint strings, that minimizes the length of one string holding other strings at constant length.

Longer term objectives:
- Define an 'exact' structure, based on the multinom library (using multinomials over integers to perform 'exact' arithmetic.) This would potentially be useful both for verifying the accuracy of multiscale solvers as well as verifying when two similar constructible points are in fact identical or slightly different.
- Implement user interfaces for basic constructions in two dimensions in SDL or SFML
- Implement user interfaces for three dimensional constructions in VR.

Contact: OpenStringGeometry@gmail.com (NOT StringGeometry@gmail.com)  
