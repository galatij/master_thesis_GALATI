function [FRIstick11,FRIstick12,FRIstick21,FRIstick22,FRIslide11,FRIslide12,FRIslide21,FRIslide22] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                    E, nu, gamma, alpha, phi, Pnu, Ptu)
    TEST = false;
    %% Extract data for face i
    n = interfData(i).normal;  % on the top face i have normal      Check: n or -n ??
    N = cpt_normal(n);
    S_n = n'*N;
    nN = repmat(n',1,8);

    t1 = interfData(i).t1;
    t2 = interfData(i).t1;
    S_t1 = t1'*N;
    S_t2 = t2'*N;
    S_t = [S_t1; S_t2];
    tT = repmat([t1, t2]', 1, 8);

    gamma = gamma/interfData(i).h;

    if (TEST)
        t = ones(3,1) - n;
        S_t = t'*N;
        tT = repmat(t', 1, 8);
        S_t = [S_t1; S_t2];
    end

    % Extract the coordinates of top and bottom faces
    loc_coo_top = coord(topol(interfData(i).etop,:),:);
    loc_coo_bot = coord(topol(interfData(i).ebottom,:),:);
    top_nod = interfData(i).top;
    bot_nod = interfData(i).bottom;

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

    % TODO: same for bottom face!!!
    ID0 = (1:3)';
    ID_top = zeros(3,1);
    ID_top(ID0~=xi_id_top) = 1:2;
    ID_bot = zeros(3,1);
    ID_bot(ID0~=xi_id_bot) = 1:2;
    
    %% Compute local contribuition on the top face (biased formulation)
    FRIstick11 = zeros(24,24);
    FRIstick12 = zeros(24,24);
    FRIstick21 = zeros(24,24);
    FRIstick22 = zeros(24,24);
    FRIslide11 = zeros(24,24);
    FRIslide12 = zeros(24,24);
    FRIslide21 = zeros(24,24);
    FRIslide22 = zeros(24,24);

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
            
            Pn = gamma*Nloc_top.*nN - (S_n*D_top*Bloc_top);                 % P 1x24
            Pnalpha = gamma*Nloc_top.*nN - alpha*(S_n*D_top*Bloc_top);
            Pt = gamma*Nloc_top.*tT - (S_t*D_top*Bloc_top);
            Ptalpha = gamma*Nloc_top.*tT - alpha*(S_t*D_top*Bloc_top);
            
            if (TEST)
                sol_test = loc_coo_top';
                sol_test = sol_test(:) + 1;
                disp(gamma*Nloc_top.*nN*sol_test);
                disp(P*sol_test);      % works correctly
            end

            % top-top (trial/test)
            FRIstick11 = FRIstick11 + 1/gamma*Ptalpha'*Pt*weights(i1)*weights(i2)*detJ;
            
            % top-bottom
            FRIstick12 = FRIstick12 + 1/gamma*Ptalpha'*(-gamma*Nloc_bot.*tT)*weights(i1)*weights(i2)*detJ;

            % bottom-top
            FRIstick21 = FRIstick21 + 1/gamma*(-gamma*Nloc_bot.*tT)'*Pt*weights(i1)*weights(i2)*detJ;

            % bottom-bottom
            FRIstick22 = FRIstick22 + 1/gamma*(-gamma*Nloc_bot.*tT)'*(-gamma*Nloc_bot.*tT)*weights(i1)*weights(i2)*detJ;
            

            %% Sliding interface

            % Prevent division by 0
            eps = 1e-10;

            Pnu_loc = Pnu(gp_idx);
            Ptu_loc = Ptu(gp_idx,:)';
            normPtu = norm(Ptu_loc);
            normPtu3 = normPtu^3 + eps;
            normPtu = normPtu + eps;

%             if (normPtu > 1e-2)
%                 ttt = Ptu(gp_idx)./normPtu^3;
%             else
%                 ttt = zeros(1,2);
%             end
%             if (normPtu > 1e-6)
%                 tt = Ptu(gp_idx)./normPtu;
%             else
%                 tt = zeros(1,2);
%             end
%             tt = tt(:);
%             ttt = ttt(:);

            wJ = weights(i1)*weights(i2)*detJ;

            trial_top = ( Pn.*(Ptu_loc/normPtu) + Pnu_loc * ( Pt/normPtu - Ptu_loc/normPtu3 * (Ptu_loc' * Pt) ) );
            
            trial_bot = ( (-gamma*Nloc_bot.*nN).*(Ptu_loc/normPtu)  ... % trial function bottom (first part)
                    + Pnu_loc * ( (-gamma*Nloc_bot.*tT)/normPtu - Ptu_loc/normPtu3 * (Ptu_loc' * (-gamma*Nloc_bot.*tT)) ) );
            
            test_top = Ptalpha;
            
            test_bot = (-gamma*Nloc_bot.*tT);

            % top-top (trial/test)
            FRIslide11 = FRIslide11 + phi/gamma*test_top'*trial_top*wJ;

            % top-bottom
            FRIslide12 = FRIslide12 + phi/gamma*test_top'*trial_bot*wJ;

            % bottom-top
            FRIslide21 = FRIslide21 + phi/gamma*test_bot'*trial_top*wJ;

            % bottom-bottom
            FRIslide22 = FRIslide22 + phi/gamma*test_bot'*trial_bot*wJ;
            
            %increment the counter
            gp_idx = gp_idx + 1;
        end
    end
    
end
