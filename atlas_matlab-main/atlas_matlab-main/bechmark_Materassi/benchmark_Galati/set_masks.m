function masksP = set_masks(Pn_gp, Pt_gp, ngauss, phi, tol_P)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Given the modified stress P evaluated at the gauss points, sets the 
%     mode (stick, slide, 0, non-smooth cases) for each gauss points.
%     The masks are used in the computation of the Jacobian and of the
%     residual.
%     masksP is a struct containing 6 entities:
%         - n0: |Pn(u)| = 0 --> first non-smooth case (triple point)
%         - nneg: Pn(u) < 0 --> [Pn(u)]_+ = 0 --> zeros Jacobian, zero residual
%         - npos: Pn(u) > 0 --> [Pn(u)]_+ = Pn(u)
%             In such a case, define Sh(u) = mu_friction*Pn(u)
%             - tstick: ||Pt(u)|| < Sh(u) --> stick
%             - t0: ||Pt(u)|| = Sh(u) --> second non-smooth case
%             - tslide: ||Pt(u)|| > Sh(u) --> project on the ball of radius Sh(u)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            normPt = norm(Pt_gp{i}(gp,:));

            if (abs(Pn_gp{i}(gp)) < tol_P)
                masksP.n0(i, gp) = true;
            elseif (Pn_gp{i}(gp) <= -tol_P)
                masksP.nneg(i, gp) = true;
            else
                masksP.npos(i, gp) = true;

                if (abs(normPt - phi * Pn_gp{i}(gp)) < tol_P)
                    masksP.t0(i, gp) = true;
                elseif ((normPt - phi * Pn_gp{i}(gp)) <= -tol_P)
                    masksP.tstick(i, gp) = true;
                else
                    masksP.tslide(i, gp) = true;
                end
            end
        end
    end
end
