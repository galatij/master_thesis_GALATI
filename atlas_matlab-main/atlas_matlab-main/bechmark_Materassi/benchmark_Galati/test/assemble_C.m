function [C] = assemble_C(ngauss,coord,ne,topol,E0,nu, ...
                          interf, interfData, gamma, alpha)
    
    if (length(E0) == 1)
        E = E0*ones(ne,1);
    else
        E = E0;
    end
    % With numerical integration of the basis functions restricted on the face

    ni = size(interf,1);
    nn = size(coord,1);
    
    k = 1;
    Clist = zeros(2*ni*144,3);
    v3 = [1;2;3];

    for i = 1 : ni
        R = interfData(i).R;
        normal = interfData(i).normal;

        % Compute contribution on the top face
        top_nod = interfData(i).top;
        D = cpt_elas_mat(E(interfData(i).etop), nu);
        Cloc_top = cpt_Cloc(ngauss, coord, topol, interfData(i).etop, top_nod, ...
            -normal, D, gamma, alpha); % 12*12
        % TODO: rotate to map on the correct side

        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        [II,JJ] = meshgrid(top_dof);
        Clist(k:k+143,:) = [JJ(:),II(:),Cloc_top(:)];
        k = k + 144;


        % Compute contribution on the bottom face
        bot_nod = interfData(i).bottom;
        D = cpt_elas_mat(E(interfData(i).ebottom), nu);
        Cloc_bottom = cpt_Cloc(ngauss, coord, topol, interfData(i).ebottom, bot_nod, ...
            normal, D, gamma, alpha); % 12*12
        % TODO: rotate to map on the correct side

        bot_dof = 3*(top_nod-1)+v3;
        bot_dof = bot_dof(:);
        [II,JJ] = meshgrid(bot_dof);
        Clist(k:k+143,:) = [JJ(:),II(:),Cloc_bottom(:)];
        k = k + 144;
        

    end
    C = sparse(Clist(:,1),Clist(:,2),Clist(:,3),3*nn,3*nn,size(Clist,1));

end
