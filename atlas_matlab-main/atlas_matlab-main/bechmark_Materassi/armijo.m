% Compute the step length with the three point parabolic model.
function [step,iarm,xp,fp,armflag] = armijo(direction,x,f0,f,maxarm,printflag)
    
    fileID = fopen("output.txt","a");
    if (fileID == -1)
        error("Cannot open the file");
    end
    sigma1 = 0.5;
    alpha = 1.e-4;
    armflag = 0;
    fp = f0;

    lambda = 1;
    lamm = 1;
    lamc = lambda;
    iarm = 0;
    step = lambda*direction;
    xt = x + step;
    ft = feval(f,xt);
    nft = norm(ft);
    nf0 = norm(f0);
    ff0 = nf0^2;
    ffc = nft^2;
    ffm = nft^2;
    while (nft >= (1 - alpha*lambda) * nf0)

        %   Apply the three point parabolic model.
        if (iarm == 0)
            lambda = sigma1*lambda;
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
        end

        % Update x; keep the books on lambda.
        step = lambda*direction;
        xt = x + step;
        lamm = lamc;
        lamc = lambda;

        % Keep the books on the function norms.
        ft = feval(f,xt);
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;

        iarm = iarm+1;
        if (iarm > maxarm)
            if (printflag)
                fprintf('*** Armijo failure, too many reductions ***\n');
                fprintf(fileID,'*** Armijo failure, too many reductions ***\n');
            end
            armflag = 1;
            xp = xt;
            fclose(fileID);
            return;
        end
        %fprintf('Armijo reduction: %e\n', lambda);
        fprintf(fileID,' ----    Armijo reduction: %e\n', lambda);
    end
    xp = xt;
    fp = ft;
    fclose(fileID);
    % end of line search
end

function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
    % Apply three-point safeguarded parabolic model for a line search.
    %
    % C. T. Kelley, April 1, 2003
    %
    % This code comes with no guarantee or warranty of any kind.
    %
    % function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
    %
    % input:
    %       lambdac = current steplength
    %       lambdam = previous steplength
    %       ff0 = value of \| F(x_c) \|^2
    %       ffc = value of \| F(x_c + \lambdac d) \|^2
    %       ffm = value of \| F(x_c + \lambdam d) \|^2
    %
    % output:
    %       lambdap = new value of lambda given parabolic model
    %
    % internal parameters:
    %       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
    %

    % Set internal parameters.
    sigma0 = 0.1;
    sigma1 = 0.5;

    % Compute coefficients of interpolation polynomial.
    % p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
    % d1 = (lambdac - lambdam)*lambdac*lambdam < 0
    %      so, if c2 > 0 we have negative curvature and default to
    %      lambdap = sigam1 * lambda.
    c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
    if (c2 >= 0)
        lambdap = sigma1*lambdac;
        return;
    end
    c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
    lambdap = -c1*0.5/c2;
    if (lambdap < sigma0*lambdac)
        lambdap = sigma0*lambdac;
    end
    if (lambdap > sigma1*lambdac)
        lambdap = sigma1*lambdac;
    end
end
