function masksP = set_masks(stress_n,stress_t, dsol, nodePairsData, gamma)

    nni = numel(nodePairsData);
    
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
            masksP.nneg(i) = true;
        elseif (abs(Pn(i)) < tol_P) % non-smooth case 1 ---> F(maskPn0) = 0, FRI(pn0,pn0) = ... 
            masksP.n0(i) = true;
        else % not open
            masksP.npos(i) = true;
            if (normPt(i) - phi*Pn(i) < -tol_P) % sticking ---> FRI(...) = sth (easy term)
                masksP.tstick(i) = true;
            elseif (abs(normPt(i) - phi*Pn(i)) < tol_P) % non-smooth case 2 --> they're almost the same, so just take the easiest one (?) 
                masksP.t0(i) = true;
            else % sliding --> FRI() = sth (difficult term)
                masksP.tslide(i) = true;
            end
        end
    end

end

%     if (TEST && norm(Pgamma_u) < 1e-12)
%         warning("Pgamma_u is zero...");
%     end