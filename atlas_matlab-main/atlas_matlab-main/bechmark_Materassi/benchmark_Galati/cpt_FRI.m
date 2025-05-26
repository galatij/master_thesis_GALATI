function [F0, FRI] = cpt_FRI(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, phi, ...
                          dsol, masksP, Pnu_gp, Ptu_gp)

    ni = size(interfData,1);
    nn = size(coord, 1);

    %% Evaluate the positive part and the projection
    v3 = [1;2;3];

    FRIstickList = zeros(ni*4*576, 3);
    FRIslideList = zeros(ni*4*576, 3);

    F0stick = zeros(3*nn, 1);
    F0slide = zeros(3*nn, 1);
    
    k = 1;
    for i = 1 : ni 

        k_base = (i-1)*4*576;

        % Compute contribution to the top face only (biased formulation)
        [FRIstick11,FRIstick12,FRIstick21,FRIstick22, ...
            FRIslide11,FRIslide12,FRIslide21,FRIslide22, ...
            RESstick_top, RESstick_bot, RESslide_top, RESslide_bot] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                        E, nu, gamma, alpha, phi, Pnu_gp{i}, Ptu_gp{i}, dsol);
        
        top_nod = topol(interfData(i).etop,:);
        bot_nod = topol(interfData(i).ebottom,:);
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);

        [II_top,JJ_top] = meshgrid(top_dof);
        FRIstickList(k_base + (1:576), :) = [JJ_top(:), II_top(:), FRIstick11(:)];

%         [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        [II_bot,JJ_top] = meshgrid(top_dof, bot_dof);
        FRIstickList(k_base + (1:576)+576, :) = [JJ_top(:), II_bot(:), FRIstick12(:)];

%         [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        [II_top,JJ_bot] = meshgrid(bot_dof, top_dof);
        FRIstickList(k_base + (1:576)+2*576, :) = [JJ_bot(:), II_top(:), FRIstick21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIstickList(k_base + (1:576)+3*576, :) = [JJ_bot(:), II_bot(:), FRIstick22(:)];


        [II_top,JJ_top] = meshgrid(top_dof);
        FRIslideList(k_base + (1:576), :) = [JJ_top(:), II_top(:), FRIslide11(:)];

%         [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        [II_bot,JJ_top] = meshgrid(top_dof, bot_dof);
        FRIslideList(k_base + (1:576)+576, :) = [JJ_top(:), II_bot(:), FRIslide12(:)];

%         [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        [II_top,JJ_bot] = meshgrid(bot_dof, top_dof);
        FRIslideList(k_base + (1:576)+2*576, :) = [JJ_bot(:), II_top(:), FRIslide21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIslideList(k_base + (1:576)+3*576, :) = [JJ_bot(:), II_bot(:), FRIslide22(:)];
        
        k = k + 576;

        % Assemble residuals for sticking
        F0stick(top_dof) = F0stick(top_dof) + RESstick_top;
        F0stick(bot_dof) = F0stick(bot_dof) + RESstick_bot;
        
        % Assemble residuals for sliding
        F0slide(top_dof) = F0slide(top_dof) + RESslide_top;
        F0slide(bot_dof) = F0slide(bot_dof) + RESslide_bot;

                
    end

    FRIstick = sparse(FRIstickList(:,2), FRIstickList(:,1), FRIstickList(:,3), 3*nn, 3*nn);   
    FRIslide = sparse(FRIslideList(:,2), FRIslideList(:,1), FRIslideList(:,3), 3*nn, 3*nn);

    [F0, FRI] = setContactMode(FRIstick, FRIslide, F0stick, F0slide, masksP, nodePairsData);

%     fprintf("||F0|| = %f\n", norm(F0));
%     fprintf("||FRI|| = %f\n", norm(FRI, 'fro'));

end

% 
%     % CHECK: ...(:,1), ...(:,2) or the reverse?
%     FRI = sparse(FRIlist11(:,2),FRIlist11(:,1),FRIlist11(:,3),3*nn,3*nn,size(FRIlist11,1)) ...
%         + sparse(FRIlist12(:,2),FRIlist12(:,1),FRIlist12(:,3),3*nn,3*nn,size(FRIlist12,1)) ...
%         + sparse(FRIlist21(:,2),FRIlist21(:,1),FRIlist21(:,3),3*nn,3*nn,size(FRIlist21,1)) ...
%         + sparse(FRIlist22(:,2),FRIlist22(:,1),FRIlist22(:,3),3*nn,3*nn,size(FRIlist22,1));
