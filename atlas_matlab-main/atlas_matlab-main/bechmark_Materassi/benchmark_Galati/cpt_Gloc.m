function [Gloc] = cpt_Gloc(ngauss, coord, topol, elem, list, ...
                                   n, gamma, D)

    TEST = false;
    loc_coo = coord(topol(elem,:),:);
    S_n = n'*cpt_normal(n);
    
    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];
    
    % Identify the fixed direction
    xi = [csi,eta,theta];
    Nloc = cpt_shape_2D(loc_coo,csi,eta,theta);
    X = ismember(Nloc*loc_coo,coord(list,:),'row');

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
    
    Gloc = zeros(24,24);
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
                sol_loc = loc_coo';
                sol_loc = sol_loc(:);
                eps = Bloc*sol_loc(:)
                stress_n = S_n*D*eps;
                tmp_coo = loc_coo';
                disp(Bloc*tmp_coo(:));        
            end
            Gloc = Gloc + (S_n*D*Bloc)'*...
                    (S_n*D*Bloc)*weights(i1)*weights(i2)*detJ;
        end
    end
    Gloc = 1/gamma*Gloc;
end
