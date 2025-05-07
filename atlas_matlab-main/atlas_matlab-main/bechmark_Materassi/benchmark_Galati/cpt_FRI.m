function [F_k, FRI] = cpt_FRI(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, phi, ...
                          dsol, stress_n, stress_t, tol_P)
    
    TEST = false;
    ni = size(interfData,1);
    nn = size(coord,1);
    nni = size(nodePairsData, 1);
    
    %% Evaluate the positive part and the projection
    v3 = [1;2;3];
    Pn = zeros(nni, 1);
    Pt = zeros(nni,2);
    % needed for the pseudo-Jacobian
    masksP = struct("n0", false(nni,1),"npos", false(nni,1),"nneg", false(nni,1),...
                   "t0", false(nni,1), "tstick", false(nni,1), "tslide", false(nni,1));

    for i = 1 : nni
        % map node -> dofs
        top_nod = nodePairsData(i).ntop;
        bot_nod = nodePairsData(i).nbottom;
        top_dof = 3*(top_nod-1)+v3;
        bot_dof = 3*(bot_nod-1)+v3;
        n = nodePairsData(i).normal;
        t1 = nodePairsData(i).t1;
        t2 = nodePairsData(i).t2;
        top_nod_loc = find([nodePairsData.ntop] == top_nod,1);
        Pn(i) = gamma*(dsol(top_dof) - dsol(bot_dof))'*n - stress_n(top_nod_loc);
        Pt(i,1) = gamma*(dsol(top_dof) - dsol(bot_dof))'*t1 - stress_t(top_nod_loc,1);
        Pt(i,2) = gamma*(dsol(top_dof) - dsol(bot_dof))'*t2 - stress_t(top_nod_loc,2);
        normPt = vecnorm(Pt,2,2);

        if (Pn(i) <= -tol_P) % open ---> F(maskPnneg) = 0
            Pn(i) = 0;
            masksP.nneg(i) = true;
        elseif (abs(Pn(i)) < tol_P) % non-smooth case 1 ---> F(maskPn0) = 0, FRI(pn0,pn0) = ... 
            Pn(i) = 0;
            masksP.n0(i) = true;
        else % not open
            masksP.npos(i) = true;
            if (normPt(i) - phi*Pn(i) < -tol_P) % sticking ---> FRI(...) = sth (easy term)
                masksP.tstick(i) = true;
            elseif (abs(normPt(i) - phi*Pn(i)) < tol_P) % non-smooth case 2 --> they're almost the same, so just take the easiest one (?) 
                masksP.t0(i) = true;
            else % sliding --> FRI() = sth (difficult term)
                masksP.tslide(i) = true;
                Pt(i,:) = Pn(i)*Pt(i,:)/normPt(i);
            end
        end
    end
    

%     if (TEST && norm(Pgamma_u) < 1e-12)
%         warning("Pgamma_u is zero...");
%     end

    FRIlist11 = zeros(ni*576, 3);
    FRIlist12 = zeros(ni*576, 3);
    FRIlist21 = zeros(ni*576, 3);
    FRIlist22 = zeros(ni*576, 3);
    k = 1;
    for i = 1 : ni 

        % Compute contribution to the top face only (biased formulation)
        [FRI11,FRI12,FRI21,FRI22] = cpt_FRIloc(ngauss, coord, topol, interfData, nodePairsData, i, ...
                    E, nu, gamma, alpha, phi, Pt, Pn, tol_P);
        
        top_nod = topol(interfData(i).etop,:);
        bot_nod = topol(interfData(i).ebottom,:);
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);

        [II_top,JJ_top] = meshgrid(top_dof);
        FRIlist11(k:k+575,:) = [JJ_top(:),II_top(:),FRI11(:)];

        [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        FRIlist12(k:k+575,:) = [JJ_top(:),II_bot(:),FRI12(:)];

        [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        FRIlist21(k:k+575,:) = [JJ_bot(:),II_top(:),FRI21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
        FRIlist22(k:k+575,:) = [II_bot(:),JJ_bot(:),FRI22(:)];
        
        k = k + 576;
                
    end
    % CHECK: ...(:,1), ...(:,2) or the reverse?
    FRI = sparse(FRIlist11(:,2),FRIlist11(:,1),FRIlist11(:,3),3*nn,3*nn,size(FRIlist11,1));
    FRI = FRI + sparse(FRIlist12(:,2),FRIlist12(:,1),FRIlist12(:,3),3*nn,3*nn,size(FRIlist12,1)) ...
              + sparse(FRIlist21(:,2),FRIlist21(:,1),FRIlist21(:,3),3*nn,3*nn,size(FRIlist21,1)) ...
              + sparse(FRIlist22(:,2),FRIlist22(:,1),FRIlist22(:,3),3*nn,3*nn,size(FRIlist22,1));
    
    F_k = FRI * dsol;
    
    % Take the positive part
    % map nni -> global dofs

    all_ntop = [nodePairsData.ntop];

    % where Pn < 0 --> set zero F, zero Jacobian
    nodes_n_neg = all_ntop(maskPnneg);
    dof_n_neg = 3*(nodes_n_neg-1)+v3;
    dof_n_neg = dof_n_neg(:);

    nodes_n_0 = all_ntop(maskPn0);
    dof_n_0 = 3*(nodes_n_0-1)+v3;
    dof_n_0 = dof_n_0(:);

    nodes_t_stick = all_ntop(maskPtstick);
    dof_t_stick = 3*(nodes_t_stick-1)+v3;
    dof_t_stick = dof_t_stick(:);

    nodes_t_0 = all_ntop(maskPt0);
    dof_t_0 = 3*(nodes_t_0-1)+v3;
    dof_t_0 = dof_t_0(:);

    % Case 1: Pn < 0
    F_k(dof_n_0 & dof_n_neg) = 0;
    FRI(dof_n_neg, dof_n_neg) = 0;

    % Case 2: Pn > 0
    % -- if ||Pt|| < f Pn --> stick --> 
    

    F_k(dof_t_neg)

    % Semi-smooth Newton where Pgamma_u = 0
    %%% e.g. take any value in the convex hull of the subdifferential, e.g.
    %%% 0.5 of the computed value
    FRI(FRIdof0, FRIdof0) = 0.5*FRI(FRIdof0, FRIdof0);
%     figure(2)
%     spy(FRI)

    if any(isnan(FRI(:))) || any(isinf(FRI(:)))
        error('FRI contains NaN or Inf values!');
    end

end

