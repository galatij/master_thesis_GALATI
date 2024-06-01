function [N,detJ] = cpt_bilinear_shape(coord,csi,eta)
%-----------------------------------------------------------------------------------------
%
% Evaluates shape 2D bilinear shape functions in a 3D space and optionally the
% Jacobian in the points defined by the 2 arrays csi and eta
%
% coord(4,3):  3D coordinates of the nodes
% csi(:):      array with first coordinate in the unit square
% eta(:):      array with second coordinate in the unit square
% N(:,4):      nodal basis functions evaluated in (csi,eta)
% detJ:        determinant of the Jacobian in (csi,eta)
%-----------------------------------------------------------------------------------------

   % Allocate room for N
   N = zeros(size(csi,1),4);
   % Evaluate shape functions
   N(:,1) = 0.25*(1.0-csi).*(1.0-eta);
   N(:,2) = 0.25*(1.0+csi).*(1.0-eta);
   N(:,3) = 0.25*(1.0+csi).*(1.0+eta);
   N(:,4) = 0.25*(1.0-csi).*(1.0+eta);

   if (nargout == 2)
      % Allocate room for derivates
      detJ = zeros(size(csi,1),1);
      dN = zeros(2,4);
      for i = 1:size(csi,1)
         % Compute derivates
         dN(1,1) = -0.25*(1 - eta(i));
         dN(2,1) = -0.25*(1 - csi(i));
         dN(1,2) =  0.25*(1 - eta(i));
         dN(2,2) = -0.25*(1 + csi(i));
         dN(1,3) =  0.25*(1 + eta(i));
         dN(2,3) =  0.25*(1 + csi(i));
         dN(1,4) = -0.25*(1 + eta(i));
         dN(2,4) =  0.25*(1 - csi(i));
         % Compute Jacobian
         V = dN*coord;
         detJ(i) = norm( cross( V(1,:) , V(2,:) ) );
      end
   end

end
