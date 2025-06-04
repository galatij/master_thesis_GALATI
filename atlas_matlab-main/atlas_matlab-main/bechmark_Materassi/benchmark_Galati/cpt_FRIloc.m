function [FRI11,FRI12,FRI21,FRI22, ...
            REStop, RESbot] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                        E, nu, gamma, alpha, phi, Pnu, Ptu, dsol, masksP)
    
    v3 = [1;2;3];

    %% Extract data for face i
    n = interfData(i).normal;   % on the top face i have normal      Check: n or -n ??
    nN = repmat(n',1,8);        % to take the normal displacement

    t1 = interfData(i).t1;
    t2 = interfData(i).t2;
    tT = repmat([t1, t2]', 1, 8);

    S_n = n' * cpt_normal(n);    % 1x6
    S_t1 = t1' * cpt_normal(n); % 1x6
    S_t2 = t2' * cpt_normal(n); % 1x6
    S_t  = [S_t1; S_t2];         % 2x6

    gamma = gamma/interfData(i).h;

    % Extract the coordinates of top and bottom faces
    loc_coo_top = coord(topol(interfData(i).etop,:),:);
    loc_coo_bot = coord(topol(interfData(i).ebottom,:),:);
    top_nod = interfData(i).top;
    bot_nod = interfData(i).bottom;

    top_dof = 3*(topol(interfData(i).etop,:)-1)+v3;
    top_dof = top_dof(:);
    bot_dof = 3*(topol(interfData(i).ebottom,:)-1)+v3;
    bot_dof = bot_dof(:);

    % Compute the elasticity tensor related to the top element
    D_top = cpt_elas_mat(E(interfData(i).etop), nu);

    %% Identify the fixed direction
    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];
    
    xi = [csi,eta,theta];
    Nloc_top = cpt_shape_2D(loc_coo_top,csi,eta,theta);
    Nloc_bot = cpt_shape_2D(loc_coo_bot,csi,eta,theta);
    X_top = ismember(Nloc_top*loc_coo_top,coord(top_nod,:),'row');
    X_bot = ismember(Nloc_bot*loc_coo_bot,coord(bot_nod,:),'row');

    for ii = 1 : 3
        if (std(xi(X_top,ii)) == 0)
            xi_id_top = ii;
            xi_val_top = mean(xi(X_top,ii));
        end
        % TODO: maybe the following is not necessary, check
        if (std(xi(X_bot,ii)) == 0)
            xi_id_bot = ii;
            xi_val_bot = mean(xi(X_bot,ii));
        end
    end

    [nodes,weights] = gausspoints(ngauss);

    ID0 = (1:3)';
    ID_top = zeros(3,1);
    ID_top(ID0~=xi_id_top) = 1:2;
    ID_bot = zeros(3,1);
    ID_bot(ID0~=xi_id_bot) = 1:2;
    
    %% Compute local contribuition on the top face (biased formulation)

    % Preallocate local matrices and res vectors
    FRI11 = zeros(24,24);
    FRI12 = zeros(24,24);
    FRI21 = zeros(24,24);
    FRI22 = zeros(24,24);
    REStop = zeros(24,1);
    RESbot = zeros(24,1);

    tmp_top = zeros(3,1);
    tmp_bot = zeros(3,1);
    tmp_top(xi_id_top) = xi_val_top;
    tmp_bot(xi_id_bot) = xi_val_bot;

    gp_idx = 1; % counter for the gauss point

    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp_top(ID_top==1) = csi;
        tmp_bot(ID_bot==1) = csi;
        for i2 = 1 : ngauss
%             fprintf("In FRIloc: i = %d, gp = %d\n", i, gp_idx);
            eta = nodes(i2);
            tmp_top(ID_top==2) = eta;
            tmp_bot(ID_bot==2) = eta;
            [Bloc_top,detJ] = cpt_shape(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top); % shape derivatives 6x24
            [Nloc_top] = cpt_shape_2D(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3)); % shape functions 1x8
            [Nloc_bot] = cpt_shape_2D(loc_coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3)); % shape functions 1x8
            Nloc_top = repelem(Nloc_top,1,3); % shape functions 1x24
            Nloc_bot = repelem(Nloc_bot,1,3);

            Pn = gamma*Nloc_top.*nN - (S_n*D_top*Bloc_top);                 % 1x24
            Pnalpha = gamma*Nloc_top.*nN - alpha*(S_n*D_top*Bloc_top);
            Pt = gamma*Nloc_top.*tT - (S_t*D_top*Bloc_top);                 % 2x24
            Ptalpha = gamma*Nloc_top.*tT - alpha*(S_t*D_top*Bloc_top);
            
            wJ = weights(i1)*weights(i2)*detJ;

            %% Sticking (RMK: F0 = FRIstick*u in case 2.1)

            if (masksP.tstick(i,gp_idx) == true ...
                    || masksP.n0(i,gp_idx) == true ...
                    || masksP.t0(i,gp_idx) == true) % stick contribution
                mode = 0;
                if (masksP.tstick(i,gp_idx) == true) % full contribution
                    mode = 1;
                elseif (masksP.n0(i,gp_idx) == true) % handle the first non-smooth case
                    mode = 1/3; % 1/3;
                elseif (masksP.t0(i,gp_idx) == true) % handle the second non-smooth case
                    mode = 0.5; % 1/2;
                end
%                 assert(mode ~= 0);

                FRI11 = FRI11 + mode*1/gamma*Ptalpha'*Pt*wJ; % top-top (test/trial)
                FRI12 = FRI12 + mode*1/gamma*Ptalpha'*(-gamma*Nloc_bot.*tT)*wJ; % top-bottom
                FRI21 = FRI21 + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*Pt*wJ; % bottom-top
                FRI22 = FRI22 + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*(-gamma*Nloc_bot.*tT)*wJ; % bottom-bottom
                
                % TODO: maybe I should compute the residual at the end
                if (masksP.tstick(i,gp_idx) == true || masksP.t0(i,gp_idx) == true)
%                     REStop = REStop + mode*1/gamma*Ptalpha'*Pt*wJ*dsol(top_dof) ...
%                         + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*Pt*wJ*dsol(bot_dof); % top
%                     RESbot = RESbot + mode*1/gamma*Ptalpha'*(-gamma*Nloc_bot.*tT)*wJ*dsol(top_dof) ...
%                         + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*(-gamma*Nloc_bot.*tT)*wJ*dsol(bot_dof); % bottom

                    REStop = REStop + mode*1/gamma*Ptalpha'*Pt*wJ*dsol(top_dof) ...
                                    + mode*1/gamma*Ptalpha'*(-gamma*Nloc_bot.*tT)*wJ*dsol(bot_dof);
                    RESbot = RESbot + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*Pt*wJ*dsol(top_dof) ...
                                    + mode*1/gamma*(-gamma*Nloc_bot.*tT)'*(-gamma*Nloc_bot.*tT)*wJ*dsol(bot_dof); % bottom
                end
            end
            if (masksP.tslide(i,gp_idx) == true ...
                    || masksP.n0(i,gp_idx) == true ...
                    || masksP.t0(i,gp_idx) == true) % sliding (RMK: F0 = FRIslide(1st part only)*u)
                mode = 0;
                if (masksP.tslide(i,gp_idx) == true) % full contribution
                    mode = 1;
                elseif (masksP.n0(i,gp_idx) == true) % handle the first non-smooth case
                    mode = 0; % 1/3
                elseif (masksP.t0(i,gp_idx) == true) % handle the second non-smooth case
                    mode = 0.5; %1/2;
                end
                
                % Prevent division by 0
                eps = 1e-10;
                
                % Evaluate Pn and Pt at the previous iter in the gauss point
                Pnu_loc = Pnu(gp_idx);
                Ptu_loc = Ptu(gp_idx,:)';
                normPtu = norm(Ptu_loc);
                normPtu3 = normPtu^3 + eps;
                normPtu = normPtu + eps;
    
                trial_top1 = Pn.*(Ptu_loc/normPtu); % trial function top (first part)
                trial_top2 = Pnu_loc * ( Pt/normPtu - Ptu_loc/normPtu3 * (Ptu_loc' * Pt) ); % trial function top (second part)
                
                trial_bot1 = (-gamma*Nloc_bot.*nN).*(Ptu_loc/normPtu); % trial function bottom (first part)
                trial_bot2 = Pnu_loc * ( (-gamma*Nloc_bot.*tT)/normPtu ...
                    - Ptu_loc/normPtu3 * (Ptu_loc' * (-gamma*Nloc_bot.*tT))  ); % trial function bottom (second part)
                
                test_top = Ptalpha;
                
                test_bot = (-gamma*Nloc_bot.*tT);
    
                % first part of the Jacobian
                FRI11 = FRI11 + mode*phi/gamma * test_top' * trial_top1 * wJ; % top-top (test/trial)
                FRI12 = FRI12 + mode*phi/gamma * test_top' * trial_top1 * wJ; % top-bottom
                FRI21 = FRI21 + mode*phi/gamma * test_bot' * trial_top1 * wJ; % bottom-top
                FRI22 = FRI22 + mode*phi/gamma * test_bot' * trial_bot1 * wJ; % bottom-bottom
                
                % residual
                if (masksP.tslide(i,gp_idx) == true || masksP.t0(i,gp_idx) == true)
                    REStop = REStop + mode*phi/gamma * test_top' * trial_top1 * wJ*dsol(top_dof) ...
                                    + mode*phi/gamma * test_top' * trial_bot1 * wJ*dsol(bot_dof); % top
                    RESbot = RESbot + mode*phi/gamma * test_bot' * trial_top1 * wJ*dsol(top_dof) ...
                                    + mode*phi/gamma * test_bot' * trial_bot1 * wJ*dsol(bot_dof); % bottom

                end
                % second part of the jacobian (accounting for dtdu)
                FRI11 = FRI11 + mode*phi/gamma*test_top'*trial_top2*wJ; % top-top (test/trial)
                FRI12 = FRI12 + mode*phi/gamma*test_top'*trial_bot2*wJ; % top-bottom
                FRI21 = FRI21 + mode*phi/gamma*test_bot'*trial_top2*wJ; % bottom-top
                FRI22 = FRI22 + mode*phi/gamma*test_bot'*trial_bot2*wJ; % bottom-bottom
            end
            % increment the counter
            gp_idx = gp_idx + 1;

        end
    end
end