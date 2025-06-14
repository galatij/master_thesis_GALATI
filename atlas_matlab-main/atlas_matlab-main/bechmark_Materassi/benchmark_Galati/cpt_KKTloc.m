function [KKT11,KKT12,KKT21,KKT22] = cpt_KKTloc(ngauss, coord, topol, interfData, i,...
    E, nu, gamma, alpha, masksP)
    
    TEST = false;
    % Extract the normal to face i
    n = interfData(i).normal;  % on the top face i have normal      Check: n or -n ??
    S_n = n'*cpt_normal(n);
    gamma = gamma/interfData(i).h;
    nN = repmat(n',1,8);

    % Extract the coordinates of top and bottom faces
    loc_coo_top = coord(topol(interfData(i).etop,:),:);
    loc_coo_bot = coord(topol(interfData(i).ebottom,:),:);
    top_nod = interfData(i).top;
    bot_nod = interfData(i).bottom;

    % Compute the elasticity tensor related to the top element
    D_top = cpt_elas_mat(E(interfData(i).etop), nu);
%     D_bot = cpt_elas_mat(E(interfData(i).etop), nu);

    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];
    
    % Identify the fixed direction
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
    
    
    % Compute local contribuition on the top face (biased formulation)
    KKT11 = zeros(24,24);
    KKT12 = zeros(24,24);
    KKT21 = zeros(24,24);
    KKT22 = zeros(24,24);
    tmp_top = zeros(3,1);
    tmp_bot = zeros(3,1);
    tmp_top(xi_id_top) = xi_val_top;
    tmp_bot(xi_id_bot) = xi_val_bot;
    gp = 1;
    
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp_top(ID_top==1) = csi;
        tmp_bot(ID_bot==1) = csi;
        for i2 = 1 : ngauss
            if (masksP.nneg(i,gp) ~= true)
                eta = nodes(i2);
                tmp_top(ID_top==2) = eta;
                tmp_bot(ID_bot==2) = eta;
                [Bloc_top,detJ] = cpt_shape(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top); % shape derivatives 6x24
                [Nloc_top] = cpt_shape_2D(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3)); % shape functions 1x8
                [Nloc_bot] = cpt_shape_2D(loc_coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3)); % shape functions 1x8
                Nloc_top = repelem(Nloc_top,1,3); % shape functions 1x24
                Nloc_bot = repelem(Nloc_bot,1,3);
                if (TEST)
    %                 disp('Shape functions (top face):');
    %                 disp(Nloc_top);
    %                 disp('Shape functions (bottom face):');
    %                 disp(Nloc_bot);
                end
                
                P = gamma*Nloc_top.*nN - (S_n*D_top*Bloc_top);                 % P 1x24
                Palpha = gamma*Nloc_top.*nN - alpha*(S_n*D_top*Bloc_top);
                if (TEST)
                    sol_test = loc_coo_top';
                    sol_test = sol_test(:) + 1;
    %                 disp(gamma*Nloc_top.*nN*sol_test);
    %                 disp(P*sol_test);      % works correctly
                    sigma = D_top*Bloc_top;
                    sigma_n = S_n*D_top*Bloc_top;
                    disp(sigma*sol_test)
                end

                % If Pn(gp) is positive, add the full contribution,
                % otherwise (Pn(gp) == 0) add half the contribution
                mode = 1;
                if (masksP.n0(i,gp) == true)
                    mode = 0.5;
                end
    
                % top-top (test/trial)
                KKT11 = KKT11 + mode*1/gamma*Palpha'*P*weights(i1)*weights(i2)*detJ;
                
                % CHECK: maybe i need to reorder the bottom dofs!!!
    
                % top-bottom
                KKT12 = KKT12 + mode*1/gamma*Palpha'*(-gamma*Nloc_bot.*nN)*weights(i1)*weights(i2)*detJ;
    
                % bottom-top
                KKT21 = KKT21 + mode*1/gamma*(-gamma*Nloc_bot.*nN)'*P*weights(i1)*weights(i2)*detJ;
    
                % bottom-bottom
                KKT22 = KKT22 + mode*1/gamma*(-gamma*Nloc_bot.*nN)'*(-gamma*Nloc_bot.*nN)*weights(i1)*weights(i2)*detJ;
            end

            gp =  gp + 1;

        end
    end
    
end
