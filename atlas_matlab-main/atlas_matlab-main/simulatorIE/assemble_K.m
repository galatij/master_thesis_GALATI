function [K,volumes] = assemble_K(ngauss,coord,topol,E,nu);
%-----------------------------------------------------------------------------------------
%
% Assembles the solid block of the stiffness matrix
% K:         stiffness matrix
% volumes:   element volumes
%
%-----------------------------------------------------------------------------------------

   nn = size(coord,1);
   ne = size(topol,1);

   % Get GAUSS point and nodes
   [nodes,weights] = gausspoints(ngauss);

   k = 1;
   ndof_ele = 3*size(topol,2);
   nnz_loc_mat = ndof_ele^2;
   Klist = zeros(ne*nnz_loc_mat,3);
   v3 = [1;2;3];
   volumes = zeros(ne,1);

   for i = 1 : ne

      % Local elastic matrix
      D = cpt_elas_mat(E(i),nu(i));

      % Init local matrix to zero
      loc_nod = topol(i,:);
      loc_coo = coord(loc_nod,:);
      Kloc = zeros(ndof_ele,ndof_ele);
      vol = 0.0;
      % Triple loop over Gauss Points
      for i1 = 1 : ngauss
         csi = nodes(i1);
         for i2 = 1 : ngauss
            eta = nodes(i2);
            for i3 = 1 : ngauss
               theta = nodes(i3);
               [Bloc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
               fac = weights(i1)*weights(i2)*weights(i3)*detJ;
               Kloc = Kloc + Bloc'*D*Bloc*fac;
               vol = vol + fac;
            end
         end
      end

      % Store volume
      volumes(i) = vol;

      % Append local entries to the global list
      loc_dof = 3*(loc_nod-1)+v3;
      loc_dof = loc_dof(:);
      [II,JJ] = meshgrid(loc_dof);
      Klist(k:k+nnz_loc_mat-1,:) = [JJ(:),II(:),Kloc(:)];
      k = k + nnz_loc_mat;

   end

   % Assemble global matrix
   K = sparse(Klist(:,1),Klist(:,2),Klist(:,3),3*nn,3*nn,size(Klist,1));

   % Force simmetry
   K = 0.5*(K + K');

end
