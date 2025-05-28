function masksP = set_masks(Pn_gp, Pt_gp, ngauss, gamma, phi, tol_P)

    ni = numel(Pn_gp);
    ngp = ngauss*2;
    v3 = [1;2;3];
    
    % needed for the pseudo-Jacobian
    masksP = struct("n0", cell(ni,1),"npos", cell(ni,1),"nneg", cell(ni,1),...
                   "t0", cell(ni,1), "tstick", cell(ni,1), "tslide", cell(ni,1));

    for i = 1 : ni
        masksP.n0{i} = false(ngp, 1);
        masksP.nneg{i} = false(ngp, 1);
        masksP.npos{i} = false(ngp, 1);
        masksP.t0{i} = false(ngp, 1);
        masksP.tstick{i} = false(ngp, 1);
        masksP.tslide{i} = false(ngp, 1);

        for gp = 1:ngp
            normPt = ...;
            if (abs(Pn_gp{i}(gp)) < tol_P) % non-smooth case 1 ---> F(maskPn0) = 0, FRI(pn0,pn0) = ... 
                masksP.n0{i}(gp) = true;
            elseif (Pn_gp{i}(gp) <= -tol_P) % open ---> F(maskPnneg) = 0
                masksP.nneg{i}(gp) = true;
            else % not open
                masksP.npos{i}(gp) = true;
                if (abs(normPt - phi*Pn_gp{i}(gp)) < tol_P) % non-smooth case 2 --> they're almost the same, so just take the easiest one (?) 
                    masksP.t0{i}(gp) = true;
                elseif (normPt - phi*Pn_gp{i}(gp) <= -tol_P) % sticking ---> FRI(...) = sth (easy term)
                    masksP.tstick{i}(gp) = true;
                else % sliding --> FRI() = sth (difficult term)
                    masksP.tslide{i}(gp) = true;
                end
            end
        end

end

%     if (TEST && norm(Pgamma_u) < 1e-12)
%         warning("Pgamma_u is zero...");
%     end