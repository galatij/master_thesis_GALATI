function [stress] = cpt_stress(ngauss,coord,topol,E,nu,sol)

    [nodes,weights] = gausspoints(ngauss);
    ne = size(topol,1);
    nn = size(coord,1);
    stress = zeros(nn, 6);
    node_count = zeros(nn,1); % counter for averaging

    v3 = [1;2;3];
    for i = 1 : ne

        % Local stiffness matrix
        D = cpt_elas_mat(E(i),nu);

        % Extract topological information
        loc_nod = topol(i,:);
        loc_coo = coord(loc_nod,:);
        loc_dof = 3*(loc_nod-1)+v3;
        loc_dof = loc_dof(:);
        loc_sol = sol(loc_dof);

        % Store stress at Gauss points (Bathe - FEProcedures, 4.3.6, p.254)
        stress_gp = zeros(ngauss^3, 6);
        gp_idx = 1;

        % Evaluate stress ath each gauss point and store in stress_gp
        for i1 = 1 : ngauss
            csi = nodes(i1);
            for i2 = 1 : ngauss
                eta = nodes(i2);
                for i3 = 1 : ngauss
                    theta = nodes(i3);
                    [Bloc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
                    % TODO: check if it's better to write the DOFs
                    % separately (loc_eps --> 6x3)
                    loc_eps = Bloc*loc_sol; % 6x24 * 24x1 --> 6x1 (or 6x3)

                    stress_gp(gp_idx,:) = (D*loc_eps)'; % 6x1 (or 6x3)
                    gp_idx = gp_idx + 1;
                end
            end
        end
        % Loop over the local nodes and add contribution to global nodes
        for n = 1 : 8
            nod = loc_nod(n);
            % TODO: take only the closest gauss points (?)
            stress(nod,:) = stress(nod,:) + mean(stress_gp,1);
            node_count(nod) = node_count(nod) + 1;
        end
    end

    for i = 1:nn
        if node_count(i) > 0
            stress(i, :) = stress(i, :) / node_count(i);
        end
    end
end
