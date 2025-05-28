function masksP = set_masks(Pn_gp, Pt_gp, ngauss, phi, tol_P)
    
    % needed for the pseudo-Jacobian
    ni = numel(Pn_gp);
    ngp = ngauss * 2;

    % Preallocate masks as logical matrices of size [ni x ngp]
    masksP.n0 = false(ni, ngp);
    masksP.nneg = false(ni, ngp);
    masksP.npos = false(ni, ngp);
    masksP.t0 = false(ni, ngp);
    masksP.tstick = false(ni, ngp);
    masksP.tslide = false(ni, ngp);

    for i = 1:ni
        for gp = 1:ngp
            normPt = norm(Pt_gp{i}(gp,:));  % Adjust indexing to match your shape

            if abs(Pn_gp{i}(gp)) < tol_P
                masksP.n0(i, gp) = true;
            elseif Pn_gp{i}(gp) <= -tol_P
                masksP.nneg(i, gp) = true;
            else
                masksP.npos(i, gp) = true;

                if abs(normPt - phi * Pn_gp{i}(gp)) < tol_P
                    masksP.t0(i, gp) = true;
                elseif normPt - phi * Pn_gp{i}(gp) <= -tol_P
                    masksP.tstick(i, gp) = true;
                else
                    masksP.tslide(i, gp) = true;
                end
            end
        end
    end
end
