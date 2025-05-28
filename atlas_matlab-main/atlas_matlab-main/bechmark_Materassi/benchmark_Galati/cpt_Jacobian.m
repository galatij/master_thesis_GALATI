function [res,J] = cpt_Jacobian(ngauss,coord,topol,E, nu,...
                                interfData, nodePairsData, K,B,alpha,gamma, ...
                                rhs,ndir,dir,state0,areaiR, ...
                                cohes,phi,tol_duT,tol_P,iter,sol)

    % Evaluate quantities at previous at gauss points k
    dsol0 = sol - state0;
    [Pn_gp, Pt_gp] = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,gamma,dsol0);          % nn*6
    masksP = set_masks(Pn_gp, Pt_gp, nodePairsData, gamma, phi, tol_P);

    % compute the residual at iteration k
    [K0, KKT] = cpt_KKT(ngauss, coord, topol, E, nu, ...
                          interfData, nodePairsData, gamma, alpha, ...
                          dsol0, masksP);

    % do the same for Friction term ...
    [F0, FRI] = cpt_FRI(ngauss, coord, topol, E, nu, ...
                         interfData, nodePairsData, gamma, alpha, phi, ...
                         dsol0, masksP, Pn_gp, Pt_gp);

    res = (K - alpha*B) * dsol0 + K0 - rhs + F0; % + F0;

    J = K - alpha*B + KKT + FRI; % + FRI;

    [J,res] = DirBC(ndir,dir,zeros(ndir,1),J,res);
    
end

