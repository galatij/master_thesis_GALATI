function [A,rhs] = DirBC(ind_dir,A,rhs);
%-----------------------------------------------------------------------------------------
%
% Applies Dirichlet conditions to Jacobian and rhs
%
%-----------------------------------------------------------------------------------------

   A(:,ind_dir) = 0;
   A = A';
   A(:,ind_dir) = 0;
   A = A';
   D = zeros(size(A,1),1);
   D(ind_dir) = 1;
   A = A + diag(sparse(D));
   rhs(ind_dir) = 0;

end
