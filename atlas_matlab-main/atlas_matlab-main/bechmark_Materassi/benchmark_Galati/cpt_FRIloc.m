function [FRIstick11,FRIstick12,FRIstick21,FRIstick22, ...
            FRIslide11,FRIslide12,FRIslide21,FRIslide22, ...
            RESstick_top, RESstick_bot, RESslide_top, RESslide_bot] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                        E, nu, gamma, alpha, phi, Pnu, Ptu, dsol)

    TEST = false;
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

    for i = 1 : 3
        if (std(xi(X_top,i)) == 0)
            xi_id_top = i;
            xi_val_top = mean(xi(X_top,i));
        end
        % TODO: maybe the following is not necessary, check
        if (std(xi(X_bot,i)) == 0)
            xi_id_bot = i;
            xi_val_bot = mean(xi(X_bot,i));
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
    FRIstick11 = zeros(24,24);
    FRIstick12 = zeros(24,24);
    FRIstick21 = zeros(24,24);
    FRIstick22 = zeros(24,24);

    FRIslide1_11 =  zeros(24,24);
    FRIslide1_12 =  zeros(24,24);
    FRIslide1_21 =  zeros(24,24);
    FRIslide1_22 =  zeros(24,24);
    
    FRIslide2_11 =  zeros(24,24);
    FRIslide2_12 =  zeros(24,24);
    FRIslide2_21 =  zeros(24,24);
    FRIslide2_22 =  zeros(24,24);

    RESstick_top = zeros(24,1);
    RESstick_bot = zeros(24,1);
    RESslide_top = zeros(24,1);
    RESslide_bot = zeros(24,1);

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
            eta = nodes(i2);
            tmp_top(ID_top==2) = eta;
            tmp_bot(ID_bot==2) = eta;
            [Bloc_top,detJ] = cpt_shape(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top); % shape derivatives 6x24
            [Nloc_top] = cpt_shape_2D(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3)); % shape functions 1x8
            [Nloc_bot] = cpt_shape_2D(loc_coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3)); % shape functions 1x8
            Nloc_top = repelem(Nloc_top,1,3); % shape functions 1x24
            Nloc_bot = repelem(Nloc_bot,1,3);

            if (TEST)
                sol_test_y = zeros(24,1);
                sol_test_z = zeros(24,1);
                for kk = 0:7
                    tmp_coo = loc_coo_top';
                    tmp_coo = tmp_coo(:);
                    sol_test_y(1+3*kk) = tmp_coo(2+3*kk);
                    sol_test_z(1+3*kk) = tmp_coo(3+3*kk);
                end
%                 sol_test = sol_test(:) + 1;
%                 disp(gamma*Nloc_top.*nN*sol_test);
%                 disp(P*sol_test);      % works correctly
                sigma = D_top*Bloc_top;
                sigma_n_y = cpt_normal(n)*D_top*Bloc_top*sol_test_y;
                sigma_t_y = S_t*D_top*Bloc_top*sol_test_y;
                sigma_n_z = cpt_normal(n)*D_top*Bloc_top*sol_test_z;
                sigma_t_z = S_t*D_top*Bloc_top*sol_test_z;
                disp(sigma*sol_test_z)
            end

            Pn = gamma*Nloc_top.*nN - (S_n*D_top*Bloc_top);                 % 1x24
            Pnalpha = gamma*Nloc_top.*nN - alpha*(S_n*D_top*Bloc_top);
            Pt = gamma*Nloc_top.*tT - (S_t*D_top*Bloc_top);                 % 2x24
            Ptalpha = gamma*Nloc_top.*tT - alpha*(S_t*D_top*Bloc_top);
            
            wJ = weights(i1)*weights(i2)*detJ;

            %% Sticking (RMK: F0 = FRIstick*u in case 2.1)

            % local Jacobian (stick)
            FRIstick11 = FRIstick11 + 1/gamma*Ptalpha'*Pt*wJ; % top-top (test/trial)
            FRIstick12 = FRIstick12 + 1/gamma*Ptalpha'*(-gamma*Nloc_bot.*tT)*wJ; % top-bottom
            FRIstick21 = FRIstick21 + 1/gamma*(-gamma*Nloc_bot.*tT)'*Pt*wJ; % bottom-top
            FRIstick22 = FRIstick22 + 1/gamma*(-gamma*Nloc_bot.*tT)'*(-gamma*Nloc_bot.*tT)*wJ; % bottom-bottom

            % residual (stick)
            RESstick_top = RESstick_top + FRIstick11*dsol(top_dof) + FRIstick21*dsol(bot_dof); % top
            RESstick_bot = RESstick_bot + FRIstick12*dsol(top_dof) + FRIstick22*dsol(bot_dof); % bottom
            

            %% Sliding (RMK: F0 = FRIslide(1st part only)*u)
            
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
            FRIslide1_11 = FRIslide1_11 + phi/gamma * test_top' * trial_top1 * wJ; % top-top (test/trial)
            FRIslide1_12 = FRIslide1_12 + phi/gamma * test_top' * trial_bot1 * wJ; % top-bottom
            FRIslide1_21 = FRIslide1_21 + phi/gamma * test_bot' * trial_top1 * wJ; % bottom-top
            FRIslide1_22 = FRIslide1_22 + phi/gamma * test_bot' * trial_bot1 * wJ; % bottom-bottom

            % residual
            RESslide_top = RESslide_top + FRIslide1_11*dsol(top_dof) + FRIslide1_21*dsol(bot_dof); % top
            RESslide_bot = RESslide_bot + FRIslide1_12*dsol(top_dof) + FRIslide1_22*dsol(bot_dof); % bottom

%             if (norm(RESslide_bot) ~= 0 || norm(RESslide_top) ~= 0)
%                 fprintf("Non-zero slide residual\n");
%             elseif (norm(RESstick_top) ~= 0 || norm(RESstick_bot) ~= 0)
%                 fprintf("Non-zero stick residual\n")
%             end
            % second part of the jacobian (accounting for dtdu)
            FRIslide2_11 = FRIslide2_11 + phi/gamma*test_top'*trial_top2*wJ; % top-top (test/trial)
            FRIslide2_12 = FRIslide2_12 + phi/gamma*test_top'*trial_bot2*wJ; % top-bottom
            FRIslide2_21 = FRIslide2_21 + phi/gamma*test_bot'*trial_top2*wJ; % bottom-top
            FRIslide2_22 = FRIslide2_22 + phi/gamma*test_bot'*trial_bot2*wJ; % bottom-bottom
            
            % increment the counter
            gp_idx = gp_idx + 1;

        end
    end
    FRIslide11 = FRIslide1_11 + FRIslide2_11;
    FRIslide12 = FRIslide1_12 + FRIslide2_12;
    FRIslide21 = FRIslide1_21 + FRIslide2_21;
    FRIslide22 = FRIslide1_22 + FRIslide2_22;
end

%             % top-top (trial/test)
%             FRIslide11 = FRIslide11 + phi/gamma*test_top'*trial_top*wJ;
% 
%             % top-bottom
%             FRIslide12 = FRIslide12 + phi/gamma*test_top'*trial_bot*wJ;
% 
%             % bottom-top
%             FRIslide21 = FRIslide21 + phi/gamma*test_bot'*trial_top*wJ;
% 
%             % bottom-bottom
%             FRIslide22 = FRIslide22 + phi/gamma*test_bot'*trial_bot*wJ;
