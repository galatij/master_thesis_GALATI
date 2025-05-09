function varargout = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,gamma,sol)
    TEST = false;

%     biased = false;

    ni = numel(interfData); % number of interface faces
    nni = numel(nodePairsData); % number of nodes at the interface
    stress = zeros(nni,6); % stress in each top node
    eps = zeros(nni,6); % strain in each top node
    node_count = zeros(nni,1); % counter for averaging

    v3 = [1;2;3];
    if nargout > 2
        Pn_gp = cell(ni,1); % for each face, store the stress in each gauss point
        Pt_gp = cell(ni,1); % for each face, interpolate the displacement (sol) in each gauss point
    end

    for i = 1:ni

        % Extract the coordinates of top and bottom faces
        % TODO: check sort(topol(....)); ???
        nod_top = topol(interfData(i).etop,:);
        nod_top_face = interfData(i).top;
        coo_top = coord(nod_top,:);
        dof_top = 3*(nod_top-1)+v3;
        dof_top = dof_top(:);
        sol_top = sol(dof_top);
        sol_top_mat = reshape(sol_top, 3, 8)';  % (8 x 3) from 24 x 1
        D_top = cpt_elas_mat(E(interfData(i).etop), nu);

        nod_bot = topol(interfData(i).ebottom,:);
        nod_bot_face = interfData(i).bottom;
        coo_bot = coord(nod_bot,:);
        dof_bot = 3*(nod_bot-1)+v3;
        dof_bot = dof_bot(:);
        sol_bot = sol(dof_bot);
        sol_bot_mat = reshape(sol_bot, 3, 8)';  % (8 x 3) from 24 x 1


        n = nodePairsData(i).normal;
        t1 = nodePairsData(i).t1;
        t2 = nodePairsData(i).t2;
        S_n = n'*cpt_normal(n);
        S_t1 = t1'*cpt_normal(t1);
        S_t2 = t2'*cpt_normal(t2);

%     if (biased)
%         D_bot = cpt_elas_mat(E(interfData(i).ebot), nu);
%     end

        % Store stress at Gauss points (Bathe - FEProcedures, 4.3.6, p.254)
        stress_gp = zeros(ngauss^2, 6);
        eps_gp = zeros(ngauss^2,6);

        Pn = zeros(ngauss^2,1);
        Pt = zeros(ngauss^2,2);

        gp_idx = 1;
    
    
        csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
        eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
        theta = [-1;-1;-1;-1;+1;+1;+1;+1];
        
        % Identify the fixed direction
        xi = [csi,eta,theta];
        Nloc_top = cpt_shape_2D(coo_top,csi,eta,theta);
        X_top = ismember(Nloc_top*coo_top,coord(nod_top_face,:),'row');
        Nloc_bot = cpt_shape_2D(coo_bot,csi,eta,theta);
        X_bot = ismember(Nloc_bot*coo_bot,coord(nod_bot_face,:),'row');

        ID0 = (1:3)';
        for j = 1 : 3
            if (std(xi(X_top,j)) == 0)
                xi_id_top = j;
                xi_val_top = mean(xi(X_top,j));
                ID_top = zeros(3,1);
                ID_top(ID0~=xi_id_top) = 1:2;
            end
            if (std(xi(X_bot,j)) == 0)
                xi_id_bot = j;
                xi_val_bot = mean(xi(X_bot,j));
                ID_bot = zeros(3,1);
                ID_bot(ID0~=xi_id_bot) = 1:2;
            end
        end

        [nodes,~] = gausspoints(ngauss);

        % Compute local contribuition on the top face (biased formulation)
        tmp_top = zeros(3,1);
        tmp_top(xi_id_top) = xi_val_top;
        tmp_bot = zeros(3,1);
        tmp_bot(xi_id_bot) = xi_val_bot;

        for i1 = 1 : ngauss
            csi = nodes(i1);
            tmp_top(ID_top==1) = csi;
            tmp_bot(ID_bot==1) = csi;
            for i2 = 1 : ngauss
                eta = nodes(i2);
                tmp_top(ID_top==2) = eta;
                tmp_bot(ID_bot==2) = eta;

                [Bloc_top,~] = cpt_shape(coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top);        % shape derivatives 6x24            [Nloc_top] = cpt_shape_2D(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3)); % shape functions 1x8
                [Nloc_top] = cpt_shape_2D(coo_top,tmp_top(1),tmp_top(2),tmp_top(3));                 % shape functions 1x8
                [Nloc_bot] = cpt_shape_2D(coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3));                 % shape functions 1x8

                loc_eps = Bloc_top * sol_top;

                stress_gp(gp_idx,:) = (D_top*loc_eps)'; % 6x1 (or 6x3)

                eps_gp(gp_idx,:) = loc_eps';

                u_gp_top = Nloc_top*sol_top_mat; % 1x8 * 8x3 --> 1x3

                u_gp_bot = Nloc_bot*sol_bot_mat; % 1x8 * 8x3 --> 1x3

                Pn(gp_idx,1) = gamma*(u_gp_top - u_gp_bot)*n - S_n*stress_gp(gp_idx,:)'; % 1x1
                Pt(gp_idx,1) = gamma*(u_gp_top - u_gp_bot)*t1 - S_t1*stress_gp(gp_idx,:)';
                Pt(gp_idx,2) = gamma*(u_gp_top - u_gp_bot)*t2 - S_t2*stress_gp(gp_idx,:)';

                gp_idx = gp_idx + 1;

                if (TEST)
                    fprintf("i = %d \t, gp: (%d, %d)\n", i, i1, i2);
                    fprintf("loc_nodes = %d %d %d %d\n", nod_top_face);
                    fprintf("loc_coords = %s\n", mat2str(coo_top(X_top,:)));
                    fprintf("sol_top = %s\n", mat2str(sol_top));
                    fprintf("loc_eps = %s\n", mat2str(loc_eps));
                end
            end
        end

        % TODO: use varargout, nargout to return stress in gauss points!
        if nargout > 2
            Pn_gp{i} = Pn;
            Pt_gp{i} = Pt;
        end

        % Loop over the local nodes and add contribution to global nodes
        for n = 1 : 4
            nod = nod_top_face(n);
            idx = find([nodePairsData.ntop] == nod, 1); 
            % TODO: take only the closest gauss points (?)
            eps(idx,:) = eps(idx,:) + mean(eps_gp,1);
            stress(idx,:) = stress(idx,:) + mean(stress_gp,1);
            node_count(idx) = node_count(idx) + 1;
        end
    end

    for n = 1:nni
        if node_count(n) > 0
            eps(n,:) = eps(n,:) / node_count(n);
            stress(n, :) = stress(n, :) / node_count(n);
        end
    end


    % Return the requested outputs
    varargout{1} = stress;
    varargout{2} = eps;
    
    if nargout > 2
        varargout{3} = Pn_gp;
        varargout{4} = Pt_gp;
    end
end
