function [F0, FRI] = cpt_FRI(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, phi, ...
                          dsol, masksP, Pnu_gp, Ptu_gp)

    ni = size(interfData,1);
    nn = size(coord, 1);

    %% Evaluate the positive part and the projection
    v3 = [1;2;3];

    FRIlist = zeros(ni*4*576, 3);
    F0 = zeros(3*nn, 1);
    
    k = 1;
    for i = 1 : ni 

        k_base = (i-1)*4*576;

        % Compute contribution to the top face only (biased formulation)
        [FRI11,FRI12,FRI21,FRI22, ...
            RES_top, RES_bot] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                        E, nu, gamma, alpha, phi, Pnu_gp{i}, Ptu_gp{i}, dsol, masksP);
        
        top_nod = topol(interfData(i).etop,:);
        bot_nod = topol(interfData(i).ebottom,:);
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);

        [II_top,JJ_top] = meshgrid(top_dof);
        FRIlist(k_base + (1:576), :) = [JJ_top(:), II_top(:), FRI11(:)];

        [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
%         [II_bot,JJ_top] = meshgrid(top_dof, bot_dof);
        FRIlist(k_base + (1:576)+576, :) = [JJ_top(:), II_bot(:), FRI12(:)];

        [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
%         [II_top,JJ_bot] = meshgrid(bot_dof, top_dof);
        FRIlist(k_base + (1:576)+2*576, :) = [JJ_bot(:), II_top(:), FRI21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIlist(k_base + (1:576)+3*576, :) = [JJ_bot(:), II_bot(:), FRI22(:)];

        k = k + 576;

        % Assemble residuals for
        F0(top_dof) = F0(top_dof) + RES_top;
        F0(bot_dof) = F0(bot_dof) + RES_bot;
        
    end

    FRI = sparse(FRIlist(:,2), FRIlist(:,1), FRIlist(:,3), 3*nn, 3*nn);   
   
%     fprintf("||F0|| = %f\n", norm(F0));
%     fprintf("||FRI|| = %f\n", norm(FRI, 'fro'));

end

% 
%     % CHECK: ...(:,1), ...(:,2) or the reverse?
%     FRI = sparse(FRIlist11(:,2),FRIlist11(:,1),FRIlist11(:,3),3*nn,3*nn,size(FRIlist11,1)) ...
%         + sparse(FRIlist12(:,2),FRIlist12(:,1),FRIlist12(:,3),3*nn,3*nn,size(FRIlist12,1)) ...
%         + sparse(FRIlist21(:,2),FRIlist21(:,1),FRIlist21(:,3),3*nn,3*nn,size(FRIlist21,1)) ...
%         + sparse(FRIlist22(:,2),FRIlist22(:,1),FRIlist22(:,3),3*nn,3*nn,size(FRIlist22,1));
