function [stress] = cpt_stress_interf(ngauss,coord,topol,interfData, nodePairsData,E,nu,sol)
    
    [nodes,weights] = gausspoints(ngauss);
    stress = zeros(nni,6);

    % Loop over the faces at interf
    for i = 1 : ni
        % Take the nodes of the top face and the normal
        loc_coo = coord(interfData(i).top);
        n = interfData(i).normal;
        S_n = n'*cpt_normal(n);
        
        csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
        eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
        theta = [-1;-1;-1;-1;+1;+1;+1;+1];
        
        % Identify the fixed direction
        xi = [csi,eta,theta];
        Nloc = cpt_shape_2D(loc_coo,csi,eta,theta);
        X = ismember(Nloc*loc_coo,coord(list,:),'row');
    
        for i = 1 : 3
            if (std(xi(X,i)) == 0)
                xi_id = i;
                xi_val = mean(xi(X,i));
            end
        end
    
        ID0 = (1:3)';
        ID = zeros(3,1);
        ID(ID0~=xi_id) = 1:2;

        D = cpt_elas_mat(E(interfData(i).etop),nu);     % 6x6 matrix
        loc_nod = interfData(i).top;                    % 4x1 list of nodes
        loc_coo = coord(loc_nod,:);
        loc_dof = 3*(loc_nod-1)+v3;
        loc_dof = loc_dof(:);
        loc_sol = sol(loc_dof);
        loc_eps = zeros(6,1);
        vol = 0.0;

        X3 = logical(kron(X', ones(1,3)));
        tmp = zeros(3,1);
        tmp(xi_id) = xi_val;

        % Computation over each face
        for i1 = 1 : ngauss
            csi = nodes(i1);
            tmp(ID==1) = csi;
            for i2 = 1 : ngauss
                eta = nodes(i2);
                tmp(ID==2) = eta;
                [Bloc,detJ] = cpt_shape(loc_coo,tmp(1),tmp(2),tmp(3),xi_id);        % shape derivatives 6x24
                loc_eps = loc_eps + D*Bloc(:,X3)*loc_sol*weights(i1)*weights(i2)*detJ;
                vol = vol + weights(i1)*weights(i2)*detJ;  
            end
        end
        loc_stress = loc_eps / vol;
        stress(i,:) = loc_stress;
    end
end