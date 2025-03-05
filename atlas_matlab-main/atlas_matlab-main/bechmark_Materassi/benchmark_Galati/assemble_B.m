function [B] = assemble_B(ngauss,coord,ne,topol,E0,nu, ...
                          interf, interfData, gamma)
    
    if (length(E0) == 1)
        E = E0*ones(ne,1);
    else
        E = E0;
    end
    % With numerical integration of the basis functions restricted on the face

    ni = size(interf,1);
    nn = size(coord,1);
    
    k = 1;
    Blist = zeros(2*ni*144,3);
    v3 = [1;2;3];

    for i = 1 : ni
        R = interfData(i).R;
        normal = interfData(i).normal;

        % Compute contribution on the top face
        top_nod = interfData(i).top;
        D = cpt_elas_mat(E(interfData(i).etop), nu);
        Gloc_top = cpt_shape_interf(ngauss, coord, topol, interfData(i).etop, top_nod, ...
            -normal, gamma, D); % 12*12
        % TODO: rotate to map on the correct side

        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        [II,JJ] = meshgrid(top_dof);
        Blist(k:k+143,:) = [JJ(:),II(:),Gloc_top(:)];
        k = k + 144;


        % Compute contribution on the bottom face
        bot_nod = interfData(i).bottom;
        D = cpt_elas_mat(E(interfData(i).ebottom), nu);
        Gloc_bottom = cpt_shape_interf(ngauss, coord, topol, interfData(i).ebottom, bot_nod, ...
            normal, gamma, D); % 12*12
        % TODO: rotate to map on the correct side

        bot_dof = 3*(top_nod-1)+v3;
        bot_dof = bot_dof(:);
        [II,JJ] = meshgrid(bot_dof);
        Blist(k:k+143,:) = [JJ(:),II(:),Gloc_bottom(:)];
        k = k + 144;

        % Map on the corresponding element

        

    end
    B = sparse(Blist(:,1),Blist(:,2),Blist(:,3),3*nn,3*nn,size(Blist,1));
    B = 0.5*(B + B');

end
