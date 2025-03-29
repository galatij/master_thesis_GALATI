function [KKTloc, KKTother] = cpt_KKTloc(ngauss, coord, topol, interfData, i, E, nu, gamma, alpha)

    % Extract the normal to face i
    n = interfData(i).normal;
    S_n = n'*cpt_normal(n);

    nN = repelem(n,4, 1);

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
    id = 0;

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
    ID = zeros(3,1);
    ID(ID0~=xi_id_top) = 1:2;
    
    % Compute local contribuition on the top face (unbiased formulation)
    KKTloc = zeros(12,12);
    KKTother = zeros(12,12);
    X3_top = logical(kron(X_top', ones(1,3)));
    X3_bot = logical(kron(X_bot', ones(1,3)));
    tmp_top = zeros(3,1);
    tmp_bot = zeros(3,1);
    tmp_top(xi_id_top) = xi_val_top;
    tmp_bot(xi_id_bot) = xi_val_bot;
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp_top(ID==1) = csi;
        tmp_bot(ID==1) = csi;
        for i2 = 1 : ngauss
            eta = nodes(i2);
            tmp_top(ID==2) = eta;
            tmp_bot(ID==2) = eta;
            [Bloc_top,detJ] = cpt_shape(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top);        % shape derivatives 6x24
            [Nloc_top] = cpt_shape_2D(loc_coo_top,tmp_top(1),tmp_top(2),tmp_top(3));                % shape functions 1x8
%             [Bloc_bot,detJ] = cpt_shape(loc_coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3),xi_id_bot);        % shape derivatives 6x24
            [Nloc_bot] = cpt_shape_2D(loc_coo_bot,tmp_bot(1),tmp_bot(2),tmp_bot(3));                % shape functions 1x8
            Nloc_top = repelem(Nloc_top(:,X_top),1,3);                                      % shape functions 1x12
            Nloc_bot = repelem(Nloc_bot(:,X_bot),1,3);
            
            
            P = Nloc_top.*nN' - gamma*(S_n*D_top*Bloc_top(:,X3_top));                                       % P 1x12
            Palpha = Nloc_top.*nN' - alpha*gamma*(S_n*D_top*Bloc_top(:,X3_top));

            % KKTloc multiplies the top displacement
            KKTloc = KKTloc + 1/gamma*P'*Palpha*weights(i1)*weights(i2)*detJ;

            % KKTother multiplies the bottom diplacement (needed for the
            % jump)
            KKTother = KKTother + 1/gamma*(-Nloc_bot.*nN')'*Palpha*weights(i1)*weights(i2)*detJ;
        end
    end
    
end
