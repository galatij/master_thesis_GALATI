function [B] = assemble_B(ngauss,coord,topol,E,nu, ...
                          interf, interfData, gamma)
    TEST = false;
    ni = size(interf,1);
    nn = size(coord,1);
    
    k = 1;
    Blist = zeros(ni*144,3);
    v3 = [1;2;3];

    for i = 1 : ni
        % Extract local data
        normal = interfData(i).normal;
        nod_top = interfData(i).top;
        etop = interfData(i).etop;
        h = interfData(i).h;
        D = cpt_elas_mat(E(etop), nu);

        % Compute contributions on the top face (biased formulation)
        [Gloc_top] = cpt_Gloc(ngauss, coord, topol, etop, nod_top, ...
            normal, gamma/h, D); % 12*12

        % Map local to global entries
        nod_elem_top = topol(etop,:);
        top_dof = 3*(nod_elem_top-1)+v3;
        top_dof = top_dof(:);
        [II,JJ] = meshgrid(top_dof);
        Blist(k:k+575,:) = [JJ(:),II(:),Gloc_top(:)];
        k = k + 576;
        
%         % Compute contribution on the bottom face (unbiased formulation)
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
