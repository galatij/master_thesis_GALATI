function [stress] = cpt_stress(ngauss,coord,elem,E,nu,sol)
% Function to compute stress on a set of elements

    ne = size(elem,1);

    [nodes,weights] = gausspoints(ngauss);

    stress = zeros(ne,6);

    v3 = [1;2;3];
    for i = 1 : ne
        % Local stiffness matrix
        D = cpt_elas_mat(E(i),nu(i));

        loc_nod = elem(i,:);
        loc_coo = coord(loc_nod,:);
        loc_dof = 3*(loc_nod-1)+v3;
        loc_dof = loc_dof(:);
        loc_sol = sol(loc_dof);
        loc_eps = zeros(6,1);
        vol = 0.0;
        for i1 = 1 : ngauss
            csi = nodes(i1);
            for i2 = 1 : ngauss
                eta = nodes(i2);
                for i3 = 1 : ngauss
                    theta = nodes(i3);
                    [Bloc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
                    loc_eps = loc_eps + Bloc*loc_sol*weights(i1)*weights(i2)*weights(i3)*detJ;
                    vol = vol + weights(i1)*weights(i2)*weights(i3)*detJ;
                end
            end
        end
        loc_stress = loc_eps / vol;
        stress(i,:) = loc_stress;
    end

end
