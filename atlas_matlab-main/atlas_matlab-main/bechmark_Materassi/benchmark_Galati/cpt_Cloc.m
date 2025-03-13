function [Cloc] = cpt_Cloc(ngauss, coord, topol, elem, list, ...
                                   n, D, gamma, alpha)

    loc_coo = coord(topol(elem,:),:);
    S_n = n'*cpt_normal(n);
    
    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];
    
    % Identify the fixed direction
    xi = [csi,eta,theta];
    Nloc = cpt_shape_2D(loc_coo,csi,eta,theta);
    X = ismember(Nloc*loc_coo,coord(list,:),'row');
    id = 0;

    for i = 1 : 3
        if (std(xi(X,i)) == 0)
            xi_id = i;
            xi_val = mean(xi(X,i));
        end
    end

    [nodes,weights] = gausspoints(ngauss);

    ID0 = [1:3]';
    ID = zeros(3,1);
    ID(ID0~=xi_id) = [1:2];
    
    % Compute local contribuition
    temp = zeros(12,12);                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cloc = zeros(12,12);
    X3 = logical(kron(X', ones(1,3)));
    tmp = zeros(3,1);
    tmp(xi_id) = xi_val;
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp(ID==1) = csi;
        for i2 = 1 : ngauss
            eta = nodes(i2);
            tmp(ID==2) = eta;
            [Bloc,detJ] = cpt_shape(loc_coo,tmp(1),tmp(2),tmp(3),xi_id);        % shape derivatives 6x24
            [Nloc] = cpt_shape_2D(loc_coo,tmp(1),tmp(2),tmp(3));                % shape functions 1x8
%             Nloc = Nloc(:,X);                                                 % shape functions 1x4 (nodes)
%             P = sum(reshape(gamma*(S_n*D*Bloc(:,X3)), 2, 3, []), 2);          % 1x12 --> 1x4 
%             Palpha= sum(reshape(alpha*gamma*(S_n*D*Bloc(:,X3)), 2, 3, []), 2); 
            Nloc = repelem(Nloc(:,X),1,3);                          % shape functions 1x12
            P = gamma*(S_n*D*Bloc(:,X3));                           % P 1x12
            Palpha = alpha*gamma*(S_n*D*Bloc(:,X3));
            Cloc = Cloc + (Nloc - P)'*(Nloc - Palpha)*weights(i1)*weights(i2)*detJ;

            [Nloc] = cpt_shape_2D(loc_coo,tmp(1),tmp(2),tmp(3));                % shape functions 1x8
            w_q = weights(i1)*weights(i2);
            Nloc = Nloc(:,X);                                                 % shape functions 1x4 (nodes)
            P = sum(reshape(gamma*(S_n*D*Bloc(:,X3)), 3, []));          % 1x12 --> 1x4 
            Palpha= sum(reshape(alpha*gamma*(S_n*D*Bloc(:,X3)), 3, []));
            % Outer product of shape functions
            N_mat = kron((Nloc - P)' * (Nloc - Palpha), eye(3));  % (3*n_nodes_per_face x 3*n_nodes_per_face)
            
            temp = temp + w_q * detJ * N_mat;
            temp - Cloc
            
        end
    end
    
end
