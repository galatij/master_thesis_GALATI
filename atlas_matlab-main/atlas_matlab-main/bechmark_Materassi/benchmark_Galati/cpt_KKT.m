function [C_k, KKT] = cpt_KKT(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, ...
                          dsol, stress_n, tol_P)
    

    ni = size(interfData,1);
    nn = size(coord,1);
    nni = size(nodePairsData, 1);
    
    %% Evaluate the modulus
    v3 = [1;2;3];
    Pgamma_u = zeros(nni, 1);
    % needed for the pseudo-Jacobian
    maskP0 = false(nni,1);
    maskPpos = false(nni,1);
    maskPneg = false(nni,1);
    for i = 1 : nni
        % map node -> dofs
        top_nod = nodePairsData(i).ntop;
        bot_nod = nodePairsData(i).nbottom;
        top_dof = 3*(top_nod-1)+v3;
        bot_dof = 3*(bot_nod-1)+v3;
        n = nodePairsData(i).normal;
        top_nod_loc = find([nodePairsData.ntop] == top_nod,1);
        Pgamma_u(i) = (dsol(top_dof) - dsol(bot_dof))'*n - stress_n(top_nod_loc);
        
        if (Pgamma_u(i) <= -tol_P)
            Pgamma_u(i) = 0;
            maskPneg(i) = true;
        elseif (abs(Pgamma_u) < tol_P)
            Pgamma_u(i) = 0;
            maskP0(i) = true;
        else
            maskPpos(i) = true;
        end
    end


    KKTlist = zeros(ni*144, 3);
    k = 1;
    for i = 1 : ni     
        % Compute contribution to the top face only (biased formulation)
        [KKTloc, KKTother] = cpt_KKTloc(ngauss, coord, topol, interfData, i, E, nu, gamma, alpha);
        
        top_nod = interfData(i).top;
        bot_nod = interfData(i).bottom;
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        [II_top,JJ_top] = meshgrid(top_dof);
        % TODO: rotate to map on the correct side
        KKTlist(k:k+143,:) = [JJ_top(:),II_top(:),KKTloc(:)];
        
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);
        [II_bot,JJ_bot] = meshgrid(top_dof, bot_dof);
        KKTlist2(k:k+143,:) = [JJ_bot(:),II_bot(:),KKTother(:)];
        k = k + 144;
        
    end
    KKT = sparse(KKTlist(:,1),KKTlist(:,2),KKTlist(:,3),3*nn,3*nn,size(KKTlist,1));
    KKT = KKT + sparse(KKTlist2(:,1),KKTlist2(:,2),KKTlist2(:,3),3*nn,3*nn,size(KKTlist2,1));

    C_k = KKT * dsol;
    
    % Take the positive part
    C_k(~maskPpos) = 0; 
    KKT(maskPneg,maskPneg) = 0;

    % Semi-smooth Newton where Pgamma_u = 0
    %%% e.g. take any value in the convex hull of the subdifferential, e.g.
    %%% 0.5 of the computed value
    KKT(maskP0,maskP0) = 0.5*KKT(maskP0, maskP0);
    
end
