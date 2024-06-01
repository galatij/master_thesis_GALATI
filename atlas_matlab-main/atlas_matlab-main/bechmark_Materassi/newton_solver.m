function [sol,iter,rnorm,conv,resvec,res] = newton_solver(f,x0,itmax,tol,maxarm,printflag,iter0)

    rnorm = 1.0;
    iter = 0;
    conv = false;
    sol = x0;
    resvec = zeros(itmax,1);
    armflag = 0;

    while (~conv)
        iter = iter + 1;

        [res,M] = f(sol,iter);

        % Solve
        x = - direct_solver(M, res);
        %@@@@@@@@@@@@@@@@@@@
        STOPPO
        %@@@@@@@@@@@@@@@@@@@

        if (iter > 1)
            [~,~,sol,~,armflag] = armijo(x,sol,res/rnorm0,@(x) f(x,iter)/rnorm0,maxarm,printflag);
        else
            sol = sol + x;
        end

        rnorm = norm(res);
        if (iter == 1)
            rnorm0 = rnorm;
        end
        resvec(iter) = rnorm/rnorm0;
        fprintf(' --- %2i %15.6e\n',iter,rnorm/rnorm0);

        conv = (rnorm < tol(1)*rnorm0+tol(2)) || iter >= itmax || armflag ~= 0;
    end
    conv = rnorm < tol(1)*rnorm0+tol(2);

    rnorm = rnorm/rnorm0;
    resvec = resvec(1:iter);
    res = res/rnorm0;

end
