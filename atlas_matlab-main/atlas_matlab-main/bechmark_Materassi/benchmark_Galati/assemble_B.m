function [B] = assemble_B(ngauss,coord,topol,E,nu, ...
                          interf, interfData, gamma)
    TEST = true;
    % With numerical integration of the basis functions restricted on the face

    ni = size(interf,1);
    nn = size(coord,1);
    
    k = 1;
    Blist = zeros(ni*144,3);
    v3 = [1;2;3];

    for i = 1 : ni
        normal = interfData(i).normal;
        
        % Compute contribution on the top face
        top_nod = interfData(i).top;
        D = cpt_elas_mat(E(interfData(i).etop), nu);
        if (TEST)
            fprintf("face: %d\ntop_nod: %d %d %d %d\n", i,top_nod);
        end
        Gloc_top = cpt_Gloc(ngauss, coord, topol, interfData(i).etop, top_nod, ...
            -normal, gamma, D); % 12*12
        % TODO: rotate to map on the correct side

        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        [II,JJ] = meshgrid(top_dof);
        Blist(k:k+143,:) = [JJ(:),II(:),Gloc_top(:)];
        k = k + 144;

%         % TODO: check if the following is needed or if to assemble only on
%         % the top side
%         % Compute contribution on the bottom face
%         bot_nod = interfData(i).bottom;
%         D = cpt_elas_mat(E(interfData(i).ebottom), nu);
%         Gloc_bottom = cpt_Gloc(ngauss, coord, topol, interfData(i).ebottom, bot_nod, ...
%             normal, gamma, D); % 12*12
%         % TODO: rotate to map on the correct side
% 
%         bot_dof = 3*(bot_nod-1)+v3;
%         bot_dof = bot_dof(:);
%         [II,JJ] = meshgrid(bot_dof);
%         Blist(k:k+143,:) = [JJ(:),II(:),Gloc_bottom(:)];
%         k = k + 144;
        

    end
    B = sparse(Blist(:,1),Blist(:,2),Blist(:,3),3*nn,3*nn,size(Blist,1));
    B = 0.5*(B + B');

end
