function [stress] = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,sol)
    ni = numel(interfData);
    nni = numel(nodePairsData);
    stress = zeros(nni,6);
    node_count = zeros(nni,1); % counter for averaging

    v3 = [1;2;3];
    
    for i = 1:ni
        % Extract the coordinates of top and bottom faces
        top_nod = interfData(i).top;
        bot_nod = interfData(i).bottom;
        coo_top = coord(topol(interfData(i).etop,:),:);
        coo_bot = coord(topol(interfData(i).ebottom,:),:);
        dof_top = 3*(top_nod-1)+v3;
        dof_top = dof_top(:);
        sol_top = sol(dof_top);

        % Store stress at Gauss points (Bathe - FEProcedures, 4.3.6, p.254)
        stress_gp = zeros(ngauss^3, 6);
        gp_idx = 1;
    
        % Compute the elasticity tensor related to the top element
        D_top = cpt_elas_mat(E(interfData(i).etop), nu);
        % D_bot = cpt_elas_mat(E(interfData(i).etop), nu);
    
        csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
        eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
        theta = [-1;-1;-1;-1;+1;+1;+1;+1];
        
        % Identify the fixed direction
        xi = [csi,eta,theta];
        Nloc_top = cpt_shape_2D(coo_top,csi,eta,theta);
        Nloc_bot = cpt_shape_2D(coo_bot,csi,eta,theta);
        X_top = ismember(Nloc_top*coo_top,coord(top_nod,:),'row');
        X_bot = ismember(Nloc_bot*coo_bot,coord(bot_nod,:),'row');

        for j = 1 : 3
            if (std(xi(X_top,j)) == 0)
                xi_id_top = j;
                xi_val_top = mean(xi(X_top,j));
            end
            % TODO: maybe the following is not necessary, check
            if (std(xi(X_bot,j)) == 0)
                xi_id_bot = j;
                xi_val_bot = mean(xi(X_bot,j));
            end
        end

        [nodes,~] = gausspoints(ngauss);
    
        ID0 = (1:3)';
        ID_top = zeros(3,1);
        ID_top(ID0~=xi_id_top) = 1:2;
        ID_bot = zeros(3,1);
        ID_bot(ID0~=xi_id_bot) = 1:2;
    
    
        % Compute local contribuition on the top face (biased formulation)
        X3_top = logical(repelem(X_top',1,3));
        X3_bot = logical(repelem(X_bot', 1,3));
        tmp_top = zeros(3,1);
        tmp_bot = zeros(3,1);
        tmp_top(xi_id_top) = xi_val_top;
        tmp_bot(xi_id_bot) = xi_val_bot;

        for i1 = 1 : ngauss
            csi = nodes(i1);
            tmp_top(ID_top==1) = csi;
            tmp_bot(ID_bot==1) = csi;
            for i2 = 1 : ngauss
                eta = nodes(i2);
                tmp_top(ID_top==2) = eta;
                tmp_bot(ID_bot==2) = eta;

                [Bloc_top,~] = cpt_shape(coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top);        % shape derivatives 6x24
                [Nloc_top] = cpt_shape_2D(coo_top,tmp_top(1),tmp_top(2),tmp_top(3));                % shape functions 1x8
                % TODO: check if it's better to write the DOFs
                % separately (loc_eps --> 6x3)
                loc_sol = sol_top.*repelem(Nloc_top(X_top),1,3)';
                loc_eps = Bloc_top(:,X3_top)*loc_sol; % 6x12 * 12x1 --> 6x1 (or 6x3)
                stress_gp(gp_idx,:) = (D_top*loc_eps)'; % 6x1 (or 6x3)
                gp_idx = gp_idx + 1;
            end
        end
        % Loop over the local nodes and add contribution to global nodes
        for n = 1 : 4
            nod = top_nod(n);
            idx = find([nodePairsData.ntop] == nod, 1); 
            % TODO: take only the closest gauss points (?)
            stress(idx,:) = stress(idx,:) + mean(stress_gp,1);
            node_count(idx) = node_count(idx) + 1;%%%%%%%%%%%%%%
        end
    end

    for n = 1:nni
        if node_count(n) > 0
            stress(n, :) = stress(n, :) / node_count(n);
        end
    end
end



% 
%     [nodes,weights] = gausspoints(ngauss);
%     ne = size(topol,1);
%     nn = size(coord,1);
%     stress = zeros(nn, 6);
%     node_count = zeros(nn,1); % counter for averaging
% 
%     v3 = [1;2;3];
%     for i = 1 : ne
% 
%         % Local stiffness matrix
%         D = cpt_elas_mat(E(i),nu);
% 
%         % Extract topological information
%         loc_nod = topol(i,:);
%         loc_coo = coord(loc_nod,:);
%         loc_dof = 3*(loc_nod-1)+v3;
%         loc_dof = loc_dof(:);
%         loc_sol = sol(loc_dof);
% 
%         % Store stress at Gauss points (Bathe - FEProcedures, 4.3.6, p.254)
%         stress_gp = zeros(ngauss^3, 6);
%         gp_idx = 1;
% 
%         % Evaluate stress ath each gauss point and store in stress_gp
%         for i1 = 1 : ngauss
%             csi = nodes(i1);
%             for i2 = 1 : ngauss
%                 eta = nodes(i2);
%                 for i3 = 1 : ngauss
%                     theta = nodes(i3);
%                     [Bloc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
%                     % TODO: check if it's better to write the DOFs
%                     % separately (loc_eps --> 6x3)
%                     loc_eps = Bloc*loc_sol; % 6x24 * 24x1 --> 6x1 (or 6x3)
% 
%                     stress_gp(gp_idx,:) = (D*loc_eps)'; % 6x1 (or 6x3)
%                     gp_idx = gp_idx + 1;
%                 end
%             end
%         end
%         % Loop over the local nodes and add contribution to global nodes
%         for n = 1 : 8
%             nod = loc_nod(n);
%             % TODO: take only the closest gauss points (?)
%             stress(nod,:) = stress(nod,:) + mean(stress_gp,1);
%             node_count(nod) = node_count(nod) + 1;
%         end
%     end
% 
%     for i = 1:nn
%         if node_count(i) > 0
%             stress(i, :) = stress(i, :) / node_count(i);
%         end
%     end
% end
