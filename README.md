# master_thesis_GALATI
Master thesis on Coulomb friction at the interface of faults.


### Modified and implemented functions:
- cpt_normal.m: returns N s.t. i can do \sigma * n as (N D Bloc), and then to get the normal stress \sigma_n i do (n' N D Bloc)
- cpt_Gloc.m and assemble_B.m: to assemble the second linear matrix
- cpt_KKTloc.m and cpt_KKT.m: to assemble the Jacobian and the residual contributions related to the KKT conditions, to be used in the Newton solver.
- cpt_FRIloc.m and cpt_FRI.m: assemble the Frictional Integral (Jacobian and Residual)
- cpt_Jacobian.m: to compute the Jacobian (and the residual), to be used inside newton_solver.m. It incorporates the linear contributions (assembled in the main) and the two non-linear terms, the first related to KKT conditions and the second related to Coulomb friction.
- cpt_stress.m: given the solution u at iteration k, compute the stress at the interface nodes (top only, biased formulation) and optionally evaluate the "modified stress" Pn, Pt at the gauss points, needed for the computation of the Jacobian. Notice that i am averaging the stress on the gauss points of the faces sharing the node, as suggested by Bathe - Finite Elements Procedures (2006).
- cpt_stress_interf.m: to take the normal traction and the tangential tractions at the interface.
- cpt_stress_tot.m: compute the stress in each element. Used only for post-processing. A refactoring may be useful since it shares mosto of the code with cpt_stress.m.
- expand_dofs.m: used in setContactMode.m. The masks for the mode (slip, stick, non-smooth case etc.) are given in each interface node. Each node is related to 3 DOFs in the matrix, but actually also other dofs contribute to the Jacobian. Hence it maps the dofs of the node to the dofs engaged in the Jacobian computation.
- setContactMode.m: after the computation of the different contributions to the Jacobian, it selects the right one (stick, slip, 0 or non-smooth case) depending on the masks computed with the previous solution. It implements the semi-smooth part of Newton.

### DONE:
1. completed the implementation of semi-smooth Newton for KKT conditions only
2. debugged cpt_normal.m, nodeManager.m, cpt_stress.m.
3. debugged the imposition of KKT condition: 
  - works properly for the mesh defined for TEST1 flag, E0 = 25000
  - does not work for the smallest mesh and for E0 = 1
  - converges in very few iterations
4. implemented the imposition of Coulomb friction

### TODO: 
1. debug FRI
2. refactory in order to modify the previous Jacobian without recomputing the whole matrices
3. store, at each time step, informations for nplas and tplas

### REMARKS:
Notice that I am using the biased formulation, for which the integral over the interface is performed only on one side of the interface. Thus I am modifying only the DOFs related to the "top (or slave)" part of the interface. The jump is accounted for by taking, through nodeManager, the value of the displacement on the two corresponding nodes at the interface, but at the end only the rows related to the TOP part are assembled. 
Another possibily, detailed by Milka (see ref.) is to use a unbiased formulation, where the same contribution is assembled on BOTH sides of the interface (multiplied by 1/2). The extension to the unbiased formulation should not be difficult once implemented the biased one, I just need to understand which implementation choice is more efficient.
