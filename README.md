# master_thesis_GALATI
Master thesis on Coulomb friction at the interface of faults.


# Modified and implemented functions:
- cpt_Gloc.m and assemble_B.m: to assemble the second linear matrix
- cpt_Cloc.m and assemble_C.m: to assemble the non-zero Jacobian, to be used in the Newton solver. At the moment i call it inside the main just for debugging on the dimensions. Since I'm getting stucked with the dimensions, I tried two possibile way to assemble the 12x12 local matrix corresponding to the local integral over the face, printing their difference to see if I got the same stuff, but it didn't work. Yet, i think that the first strategy (the same that I used in cpt_Gloc) may be the correct one, if any.
- cpt_normal.m: returns N s.t. i can do \sigma * n as (N D Bloc), and then to get the normal stress \sigma_n i do (n' N D Bloc)
