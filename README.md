# master_thesis_GALATI
Master thesis on Coulomb friction at the interface of faults.


### Modified and implemented functions:
- cpt_normal.m: returns N s.t. i can do \sigma * n as (N D Bloc), and then to get the normal stress \sigma_n i do (n' N D Bloc)
- cpt_Gloc.m and assemble_B.m: to assemble the second linear matrix
- cpt_KKTloc.m and cpt_KKT.m: to assemble the Jacobian and the residual contributions related to the KKT conditions, to be used in the Newton solver. TODO:
  - choose where to cpt sigman(u) and where to store it (if needed) for the computation of Pgamma(u);
  - modify KKT matrix (for the Jacobian) depending of maskP, applying a choice for semi-smooth Newton
- cpt_Jacobian.m: to compute the Jacobian (and the residual), to be used inside newton_solver.m. It will incorporate the linear contributions (assembled in the main) and the two non-linear terms, the first related to KKT conditions and the second related to Coulomb friction.
- cpt_stress_interf.m: to compute the stress at the interface, given the solution u at iteration k. Needed for the computation of Pgamma(u) in both the nonlinear terms.ù

### TODO: 
1. choose where to compute/store sigma(u) at the interface
2. complete the implementation of semi-smooth Newton for KKT conditions only
3. debug and test this easier problem
4. keep going with the friction term (it should be similar to the KKT term)

### REMARKS:
Notice that I am using the unbiased formulation, for which the integral over the interface is performed only on one side of the interface. Thus I am modifying only the DOFs related to the "top (or slave)" part of the interface. The jump is accounted for by taking, through nodeManager, the value of the displacement on the two corresponding nodes at the interface, but at the end only the rows related to the TOP part are assembled. 
Another possibily, detailed by Milka (see ref.) is to use a unbiased formulation, where the same contribution is assembled on BOTH sides of the interface (multiplied by 1/2).
