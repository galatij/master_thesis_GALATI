function interf2e = map_IE_elem(interf,elem)
%-----------------------------------------------------------------------------------------
%
% Function to create a mapping between IEs and finite elements:
% interf2e is a sparse matrix where the row index corresponds to the IE number and the
% column index corresponds to the adjacent finite element. The entry is 1 for the BOTTOM
% and 2 for the TOP
%
% NOTE: this procedure is really inefficient
%
%-----------------------------------------------------------------------------------------

   ni = size(interf,1);
   ne = size(elem,1);
   face = [1,2,3,4;1,2,6,5;2,3,7,6;4,3,7,8;1,4,8,5;5,6,7,8];
   irow = zeros(2*ni,1);
   jcol = zeros(2*ni,1);
   coef = zeros(2*ni,1);
   intB = sort(interf(:,1:4),2);
   intT = sort(interf(:,5:8),2);
   k = 0;
   for i = 1:ne
       for j = 1:6
           f = sort(elem(i,face(j,:)));
           [flag, ii] = ismember(f, intB, 'rows');
           if (flag)
               k = k + 1;
               irow(k) = ii;
               jcol(k) = i;
               coef(k) = 1;
           else
              [flag, ii] = ismember(f, intT, 'rows');
              if (flag)
                  k = k + 1;
                  irow(k) = ii;
                  jcol(k) = i;
                  coef(k) = 2;
              end
           end
       end
   end
   interf2e = sparse(irow,jcol,coef);

end
