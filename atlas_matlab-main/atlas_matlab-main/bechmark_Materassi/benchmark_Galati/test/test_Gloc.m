function [sigma_n] = test_Gloc(ngauss, coord, topol, elem, list, ...
                                   n, D)

    loc_coo = coord(topol(elem,:),:);
    S_n = n'*cpt_normal(n);
    
    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];
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
    ID0 = [1:3]';
    ID = zeros(3,1);
    ID(ID0~=xi_id) = [1:2];    
    X3 = logical(kron(X', ones(1,3)));
    tmp = zeros(3,1);
    tmp(xi_id) = xi_val;

    [nodes,~] = gausspoints(ngauss);
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp(ID==1) = csi;
        for i2 = 1 : ngauss
            eta = nodes(i2);
            tmp(ID==2) = eta;
            [Bloc,~] = cpt_shape(loc_coo,tmp(1),tmp(2),tmp(3),xi_id);
            B_u = repmat([diag(ones(3,1)); zeros(3,3)],1,4);
            sigma_n = sigma_n + (S_n*D*Bloc(:,X3))';
        end
    end
end
