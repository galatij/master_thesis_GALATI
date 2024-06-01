function C = cpt_global_stab(K, Bt);
%-----------------------------------------------------------------------------------------
%
% Function to compute the "global stabilization" block from K and Bt
% CHIEDERE AD ANDREA MEGLIO
%
%-----------------------------------------------------------------------------------------

   % Retrieve number of Lagrange multipliers and init an empty C block
   nl = size(Bt,1);
   C = sparse(nl,nl);

   % Count how many nodes IE share
   Bt_nodal = spones(Bt(1:3:end,1:3:end));
   BtB_nodal = Bt_nodal*Bt_nodal';

   I3 = eye(3);
   Z3 = zeros(3);
   signMat = [I3,Z3;Z3,-I3];
   % Loop over the interface elements
   for i = 1:size(BtB_nodal,1)
      % Get the list of IE elements sharing nodes with i
      [ind,~,bb] = find(BtB_nodal(:,i));
      if (length(ind) > 0)
         % Get the diagonal element (max number of connections)
         diagBtB = bb(ind==i);
         % List of IE sharing one edge
         faceLoc = ind(bb==diagBtB/2);
         % Trick to manage boundaries
         if (diagBtB ~= max(bb) || length(faceLoc) == 0)
            faceLoc = ind(bb >= diagBtB);
         end
         % To avoid to double assemble the same entries
         faceLoc = faceLoc(faceLoc>i);
         dof1 = (3*(i-1)+1) : (3*i);
         [~,nodes1,~] = find(Bt(dof1,:));
         for j = 1:length(faceLoc)
            dof2 = (3*(faceLoc(j)-1)+1) : (3*faceLoc(j));
            [~, nodes2, ~] = find(Bt(dof2,:));
            dofLoc = intersect(nodes1, nodes2, 'sorted');
            BtSign = signMat*Bt([dof2,dof1],dofLoc);
            CLoc = BtSign * diag(diag(K(dofLoc,dofLoc)))^-1 * BtSign';
            C([dof1,dof2],[dof1,dof2]) = C([dof1,dof2],[dof1,dof2]) + CLoc;
         end
      end
   end

   % Ensure symmetry
   C = 0.5*(C+C');

end
