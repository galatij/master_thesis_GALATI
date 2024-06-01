function C = cpt_global_stab(K, Bt);

    ni = size(Bt,1);
    C = sparse(ni,ni);

    % Jump penalization
    B1t = spones(Bt);
    BtB = B1t * B1t';
    BtB_nod = compress_mat(BtB, 3);

    I3 = eye(3);
    Z3 = zeros(3);
    signMat = [I3,Z3;Z3,-I3];
    for i = 1 : size(BtB_nod,1)
        [ind, ~, bb] = find(BtB_nod(:,i));
        if (length(ind) > 0)
            diagBtB = bb(ind==i);
            faceLoc = ind(bb==diagBtB/2);
            % Trick to manage boundaries
            if (diagBtB ~= max(bb) || length(faceLoc) == 0)
                faceLoc = ind(bb>=diagBtB);
            end
            % To avoid to double assemble the same entries
            faceLoc = faceLoc(faceLoc>i);
            dof1 = (3*(i-1)+1) : (3*i);
            [~, nodes1, ~] = find(Bt(dof1,:));
            for j = 1 : length(faceLoc)
                dof2 = (3*(faceLoc(j)-1)+1) : (3*faceLoc(j));
                [~, nodes2, ~] = find(Bt(dof2,:));
                dofLoc = intersect(nodes1, nodes2, 'sorted');
                BtSign = signMat * Bt([dof2,dof1],dofLoc);
                CLoc = BtSign * diag(diag(K(dofLoc,dofLoc)))^-1 * BtSign';
                C([dof1,dof2],[dof1,dof2]) = C([dof1,dof2],[dof1,dof2]) + CLoc;
            end
        end
    end

    C = 0.5*(C+C');
end
