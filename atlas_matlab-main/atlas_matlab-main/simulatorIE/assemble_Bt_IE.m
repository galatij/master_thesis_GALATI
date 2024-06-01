function [Bt] = assemble_Bt_IE(ngauss,coord,interf,interfData)
%-----------------------------------------------------------------------------------------
%
% Function to compute and assemble the contribution of the Lagrange multiplier to the
% stiffness matrix (Bt block) with numerical integration of the basis functions on the
% IE face
% ngauss:         number of gauss points
% coord(:,3):     global coordinates of the mesh
% interf(:,8):    nodal connections of the IE element
% interfData(:)   structure collecting IE information
%
%-----------------------------------------------------------------------------------------

   % Retrieve Gauss points and weights
   [Gnodes,Gweights] = gausspoints(ngauss);

   ni = size(interf,1);
   nn = size(coord,1);
   ndof_ie = 3*size(interf,2);
   nlambda_ie = 3;
   nnz_loc_mat = ndof_ie*nlambda_ie;
   irow = zeros(nnz_loc_mat*ni,1);
   jcol = zeros(nnz_loc_mat*ni,1);
   coef = zeros(nnz_loc_mat*ni,1);
   v3 = [1;2;3];
   k = 1;
   for i = 1:ni

      % Get rotation matrix
      R_l2g = interfData(i).R_l2g;
      % Get nodal coordinates of the face
      nod_coor = coord(interf(i,1:4),:);
      % Compute nodal area of the face
      int_area = cpt_area_int(Gnodes,Gweights,nod_coor);

      % Assemble contribution of the top
      list = interfData(i).top;
      for j = 1 : 4
         irow(k:k+2) = 3*(i-1)+1;
         jcol(k:k+2) = 3*(list(j)-1)+v3;
         irow(k+3:k+5) = 3*(i-1)+2;
         jcol(k+3:k+5) = 3*(list(j)-1)+v3;
         irow(k+6:k+8) = 3*(i-1)+3;
         jcol(k+6:k+8) = 3*(list(j)-1)+v3;
         coef(k:k+8) = int_area(j)*R_l2g(:);
         k = k + 9;
      end

      % Assemble contribution of the bottom
      list = interfData(i).bottom;
      for j = 1 : 4
         irow(k:k+2) = 3*(i-1)+1;
         jcol(k:k+2) = 3*(list(j)-1)+v3;
         irow(k+3:k+5) = 3*(i-1)+2;
         jcol(k+3:k+5) = 3*(list(j)-1)+v3;
         irow(k+6:k+8) = 3*(i-1)+3;
         jcol(k+6:k+8) = 3*(list(j)-1)+v3;
         coef(k:k+8) = -int_area(j)*R_l2g(:);
         k = k + 9;
      end
   end
   Bt = sparse(irow,jcol,coef,3*ni,3*nn);

end
