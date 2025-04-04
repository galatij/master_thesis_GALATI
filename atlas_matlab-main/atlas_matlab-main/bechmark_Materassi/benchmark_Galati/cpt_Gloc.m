function [Gloc] = cpt_Gloc(ngauss, coord, topol, elem, list, ...
                                   n, gamma, D)

    TEST = false;
    % list has the 4 global indices defining the face
    % coordinates of the top element (cell) (8x3 matrix)
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

    % Loop over the 3 ref. coords to determine the constant one
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
    if (TEST)
        Gloc_cmp = zeros(12,12);
    end
    Gloc = zeros(12,12);
    X3 = logical(kron(X', ones(1,3)));     % X = repmat(X',1,3);   TODO: check storage of components
%     e.g. on command line
%     X = [1;0;1; 1];
%     X3 = logical(kron(X',ones(1,3)))  % --> [1,1,1,0,0,0,1,1,1,1,1,1]
    tmp = zeros(3,1);
    tmp(xi_id) = xi_val;
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp(find(ID==1)) = csi;
        for i2 = 1 : ngauss
            eta = nodes(i2);
            tmp(find(ID==2)) = eta;
            [Bloc,detJ] = cpt_shape(loc_coo,tmp(1),tmp(2),tmp(3),xi_id);
            if (TEST)
%                 Bloc(:,X)
%                 D*Bloc(:,X)
%                 S_n*D*Bloc(:,X)
                B_u = repmat([diag(ones(3,1)); zeros(3,3)],1,4);
                Gloc_cmp = Gloc_cmp + (S_n*D*B_u)'*...
                            (S_n*D*Bloc(:,X3))*weights(i1)*weights(i2)*detJ;
%                 [Btmp,detJ] = cpt_shape(loc_coo,0,0,0,xi_id);
%                 tmp_coo = loc_coo(X,:);
%                 Btmp = Btmp(:,X3);
%                 Btmp*tmp_coo(:);

        
            end
            Gloc = Gloc + (S_n*D*Bloc(:,X3))'*...
                    (S_n*D*Bloc(:,X3))*weights(i1)*weights(i2)*detJ;
        end
    end
    Gloc = 1/gamma*Gloc;

    if (TEST)
        Gloc_cmp = 1/gamma*Gloc_cmp;
    end

end
