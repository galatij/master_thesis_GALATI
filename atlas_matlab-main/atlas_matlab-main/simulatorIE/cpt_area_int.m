function [int_area] = cpt_area_int(Gnodes,Gwgts,loc_coo)
%-----------------------------------------------------------------------------------------
%
% Function to compute nodal area for a input face
% Gnodes:        list of 1D gauss coordinates
% Gwgts:         list of 1D gauss weights
% loc_coo(:,3):  3D coordinates of the nodes (size(loc_coo,1) gives the number of nodes)  
% int_area(:):   area associated to each node
%
%-----------------------------------------------------------------------------------------

   n_gauss = numel(Gnodes);

   % Evaluate 2D integration points and weights
   csi = kron(ones(n_gauss,1),Gnodes);
   eta = kron(ones(1,n_gauss),-Gnodes');
   area_wgt = Gwgts'*Gwgts;
   % Store integration points and weights in 1D array
   csi = csi(:);
   eta = eta(:);
   area_wgt = area_wgt(:);

   % Compute basis functions and Jacobian
   [N_nodal,detJ] = cpt_bilinear_shape(loc_coo,csi,eta);

   % Compute nodal area
   int_area = (detJ.*area_wgt)'*N_nodal;

end
