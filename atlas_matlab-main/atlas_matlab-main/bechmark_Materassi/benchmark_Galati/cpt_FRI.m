function [F_k, FRI] = cpt_FRI(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, phi, ...
                          dsol, masksP, Pnu_gp, Ptu_gp)

    ni = size(interfData,1);

    %% Evaluate the positive part and the projection
    v3 = [1;2;3];


    FRIstickList11 = zeros(ni*576, 3);
    FRIstickList12 = zeros(ni*576, 3);
    FRIstickList21 = zeros(ni*576, 3);
    FRIstickList22 = zeros(ni*576, 3);

    FRIslideList11 = zeros(ni*576, 3);
    FRIslideList12 = zeros(ni*576, 3);
    FRIslideList21 = zeros(ni*576, 3);
    FRIslideList22 = zeros(ni*576, 3);
    
    k = 1;
    for i = 1 : ni 

        % Compute contribution to the top face only (biased formulation)
        [FRIstick11,FRIstick12,FRIstick21,FRIstick22, ...
            FRIslide11,FRIslide12,FRIslide21,FRIslide22] = ...
                    cpt_FRIloc(ngauss, coord, topol, interfData, i, ...
                        E, nu, gamma, alpha, phi, Pnu_gp{i}, Ptu_gp{i});
        
        top_nod = topol(interfData(i).etop,:);
        bot_nod = topol(interfData(i).ebottom,:);
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);

        [II_top,JJ_top] = meshgrid(top_dof);
        FRIstickList11(k:k+575,:) = [JJ_top(:),II_top(:),FRIstick11(:)];

        [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        FRIstickList12(k:k+575,:) = [JJ_top(:),II_bot(:),FRIstick12(:)];

        [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        FRIstickList21(k:k+575,:) = [JJ_bot(:),II_top(:),FRIstick21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIstickList22(k:k+575,:) = [II_bot(:),JJ_bot(:),FRIstick22(:)];


        [II_top,JJ_top] = meshgrid(top_dof);
        FRIslideList11(k:k+575,:) = [JJ_top(:),II_top(:),FRIslide11(:)];

        [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        FRIslideList12(k:k+575,:) = [JJ_top(:),II_bot(:),FRIslide12(:)];

        [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        FRIslideList21(k:k+575,:) = [JJ_bot(:),II_top(:),FRIslide21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIslideList22(k:k+575,:) = [II_bot(:),JJ_bot(:),FRIslide22(:)];
        
        k = k + 576;
                
    end

    % Combine lists
    FRIstickList = [FRIstickList11; FRIstickList12; FRIstickList21; FRIstickList22];
    FRIslideList = [FRIslideList11; FRIslideList12; FRIslideList21; FRIslideList22];
    
    % Create masks (4 blocks per interface)
    mask_stick = repmat(masksP.tstick, 576 * 4, 1);
    mask_slide = repmat(masksP.tslide, 576 * 4, 1);
    
    % Apply masks
    FRIstickList(~mask_stick,:) = 0;
    FRIslideList(~mask_slide,:) = 0;
    
    % Sum and clean
    FRIall = FRIstickList + FRIslideList;
    FRIall = FRIall(FRIall(:,3) ~= 0, :);
    
    % Assemble sparse matrix
    ndof = 3 * size(coord,1);
    F_k = sparse(FRIall(:,2), FRIall(:,1), FRIall(:,3), ndof, ndof);

    

end

% 
%     % CHECK: ...(:,1), ...(:,2) or the reverse?
%     FRI = sparse(FRIlist11(:,2),FRIlist11(:,1),FRIlist11(:,3),3*nn,3*nn,size(FRIlist11,1)) ...
%         + sparse(FRIlist12(:,2),FRIlist12(:,1),FRIlist12(:,3),3*nn,3*nn,size(FRIlist12,1)) ...
%         + sparse(FRIlist21(:,2),FRIlist21(:,1),FRIlist21(:,3),3*nn,3*nn,size(FRIlist21,1)) ...
%         + sparse(FRIlist22(:,2),FRIlist22(:,1),FRIlist22(:,3),3*nn,3*nn,size(FRIlist22,1));
%     
%     F_k = FRI * dsol;
%     
%     % Take the positive part
%     % map nni -> global dofs
% 
%     all_ntop = [nodePairsData.ntop];
% 
%     % where Pn < 0 --> set zero F, zero Jacobian
%     nodes_n_neg = all_ntop(masksP.nneg);
%     dof_n_neg = 3*(nodes_n_neg-1)+v3;
%     dof_n_neg = dof_n_neg(:);
% 
%     nodes_n_0 = all_ntop(masksP.n0);
%     dof_n_0 = 3*(nodes_n_0-1)+v3;
%     dof_n_0 = dof_n_0(:);
% 
%     nodes_t_stick = all_ntop(masksP.stick);
%     dof_t_stick = 3*(nodes_t_stick-1)+v3;
%     dof_t_stick = dof_t_stick(:);
% 
%     nodes_t_0 = all_ntop(masksP.t0);
%     dof_t_0 = 3*(nodes_t_0-1)+v3;
%     dof_t_0 = dof_t_0(:);
% 
%     % Case 1: Pn < 0
%     F_k(dof_n_0 & dof_n_neg) = 0;
%     FRI(dof_n_neg, dof_n_neg) = 0;
% 
%     % Case 2: Pn > 0
%     % -- if ||Pt|| < f Pn --> stick --> 
%     
% 
%     F_k(dof_t_neg)
% 
%     % Semi-smooth Newton where Pgamma_u = 0
%     %%% e.g. take any value in the convex hull of the subdifferential, e.g.
%     %%% 0.5 of the computed value
%     FRI(FRIdof0, FRIdof0) = 0.5*FRI(FRIdof0, FRIdof0);
% %     figure(2)
% %     spy(FRI)
% 
%     if any(isnan(FRI(:))) || any(isinf(FRI(:)))
%         error('FRI contains NaN or Inf values!');
%     end
% 
% end

