function [sol, convNRvec, convNR, ittot] = ...
    solve_NL_CM(K, B, alpha, gamma, E, nu, rhs, interfData, nodePairsData, areaiR, sol0, state0,...
                ndir, dir, dirval, noConvItMax, itmax_NR, tol_NR, ...
                cohes, phi, tol_sig, tol_duNc, tol_duT, tol_P, ngauss, coord, topol, ...
                volumes)

    maxarm = 3;
    printflag = true;

    sol = sol0;     % previous time step
    sol(dir(1:ndir)) = dirval(1:ndir);
    convNR = false;
    exitFlag = false;
    ittot = 0;
    itNoConv = 0;
    convNRvec = zeros(itmax_NR,1);

    while (~exitFlag)

        f = @(sol,iter) cpt_Jacobian(ngauss,coord,topol,E,nu,...
                                interfData,nodePairsData,K,B,alpha,gamma, ...
                                rhs,ndir,dir,state0,areaiR, ...
                                cohes,phi,tol_duT,tol_P,iter,sol);

        % TODO: choose where to compute the stress
        % notice that the projection of the stress are computed in the
        % integrals of Nitsche method

        [solNew,iter,rnorm,convNR,resvec] = newton_solver(f,sol,itmax_NR,tol_NR,maxarm,printflag);
		

        ittot = ittot + iter;
        convNRvec(ittot-iter+1:ittot) = resvec;

        if (convNR)
            sol = solNew;
            exitFlag = true;
            fprintf("Newton solver has converged with %d iterations\n", iter);
            % TODO: check nplas tplas and tplasnew (just for output)
        else
            itNoConv = itNoConv + 1;
        end

        exitFlag = exitFlag || (itNoConv > noConvItMax);
    end

    convNRvec = convNRvec(1:ittot);

end
