function [sol_u, sol_l, nplas, tplas, convNRvec, convFlag, C, ittot] = ...
    solve_NL_CM(nn, ni, K, Bt, rhs, interfData, areaiR, sol0, state0, nplas, tplas, ...
                ndir, dir, dirval, actSetItMax, noConvItMax, itmax_NR, tol_NR, ...
                cohes, phi, tol_sig, tol_duNc, tol_duT, iter0, ngauss, coord, topol, ...
                volumes, edgeData, f2e, E, nu)
    fileID = fopen("output.txt","a");
    if (fileID == -1)
        error("Cannot open the file");
    end
    maxarm = 3;
    printflag = true;

    indU = (1:3*nn);
    indL = (3*nn) + (1:3*ni);

    B = Bt';

    Corig = cpt_global_stab(K, Bt);
    Corig = Corig - diag(sum(Corig,1));

    C = setStabilizationMat(Corig, nplas, tplas);

    sol = sol0;
    sol(dir(1:ndir)) = dirval(1:ndir);
    convNR = false;
    exitFlag = false;
    checkActiveSet = false;
    ittot = 0;
    itNoConv = 0;
    convNRvec = zeros(actSetItMax*itmax_NR,1);
    trialStep = 0;

    % If there are no slipping elements, it's already an elasticStep!
    if (sum(tplas(~nplas)) == 0)
        trialStep = 1;
    end

    % Active set variables
    itA = 0;
    statusL = zeros(2*ni,itmax_NR);
    skipActiveSet = false;
    tplasnew = false(ni,1);

    while (~exitFlag)

        f = @(x,iter) cpt_residual(K,B,Bt,C,rhs,ndir,dir,state0,sol0,indU,indL,areaiR, ...
                                   nplas,tplas,tplasnew,cohes,phi,tol_duT,interfData,iter,x);

        [solNew,iter,rnorm,convNR,resvec] = newton_solver(f,sol,itmax_NR,tol_NR,maxarm,printflag,ittot);
		STOPPO
        ittot = ittot + iter;
        convNRvec(ittot-iter+1:ittot) = resvec;

        fprintf('%6i %15.6e (%5i %5i)\n', iter0+ittot, rnorm, sum(nplas), sum(tplas));
        fprintf(fileID, '%6i %15.6e (%5i %5i)\n', iter0+ittot, rnorm, sum(nplas), sum(tplas));

        if (~skipActiveSet)
            [nplas,tplas,tplasnew,checkActiveSet,skipActiveSet,itA,statusL] = ...
                activeSet(Bt,indU,indL,nplas,tplas,itA,statusL,areaiR,cohes,phi,tol_duNc, ...
                          tol_sig,convNR,solNew);
            C = setStabilizationMat(Corig, nplas, tplas);
        end

        if (convNR)
            sol = solNew;
            exitFlag = checkActiveSet;
        elseif (skipActiveSet)
            sol = solNew;
            exitFlag = true;
        elseif (trialStep == 0)
            % elasticStep
            exitFlag = false;
            tplas = nplas;
            trialStep = trialStep + 1;
            fprintf(' ------------- elastic step -------------\n');
            fprintf(fileID, ' ------------- elastic step -------------\n');
            C = setStabilizationMat(Corig, nplas, tplas);
        else
            exitFlag = true;
        end

        if (~convNR)
            itNoConv = itNoConv + 1;
        end

        exitFlag = exitFlag || (itA > actSetItMax) || (itNoConv > noConvItMax);
    end

    convFlag = convNR && checkActiveSet && (itA <= actSetItMax);

    sol_u = sol(indU);
    sol_l = sol(indL);
    convNRvec = convNRvec(1:ittot);

    % --------------------------------------------------------------
    % JUST FOR OUTPUT
    % --------------------------------------------------------------
    C = setStabilizationMat(Corig, nplas, tplas);
    % --------------------------------------------------------------
    % JUST FOR OUTPUT
    % --------------------------------------------------------------
    fclose(fileID);
end

function [res,J] = cpt_residual(K,B,Bt,C,rhs,ndir,dir,state0,sol0,indU,indL,areaiR, ...
                                nplas,tplas,tplasnew,cohes,phi,tol_duT,interfData,iter,sol)

    cptJ = false;
    if (nargout > 1)
        cptJ = true;
    end

    I2 = speye(2);
    nn = length(indU)/3;
    ni = length(indL)/3;

    if (cptJ)
        % Remove constraints for open and sliding elements
        BtDir = setCouplingMat(Bt, nplas, tplas);
    end

    % Previous sliding
    IDn = 3*((1:ni)-1)+1;
    dsol0 = sol0 - state0;
    gT0 = Bt*dsol0(indU);
    gT0(IDn) = 0;

    Cl0 = C*dsol0(indL);
    Cl0(IDn) = 0;

    dsol = sol - state0;
    res_u = (K*dsol(indU) + B*dsol(indL)) - rhs(indU);
    res_l = (Bt*dsol(indU) - C*dsol(indL) - (gT0 - Cl0)) - rhs(indL);

    if (cptJ)
        indE = 0;
        irowE = zeros(48*ni,1);
        jcolE = zeros(48*ni,1);
        valsE = zeros(48*ni,1);

        indF = 0;
        irowF = zeros(3*ni,1);
        jcolF = zeros(3*ni,1);
        valsF = zeros(3*ni,1);
    end

    % Relative displacements (w/o stabilization effects)
    dulocTot0 = areaiR*(Bt*dsol0(indU));
    dulocTot = areaiR*(Bt*dsol(indU));
    sol_l = sol(indL);

    for i = 1 : ni
        if (nplas(i))
            idof = 3*(i-1)+(1:3);
            res_l(idof(1)) = interfData(i).area * sol_l(idof(1));
            res_l(idof(2)) = interfData(i).area * sol_l(idof(2));
            res_l(idof(3)) = interfData(i).area * sol_l(idof(3));
            if (cptJ)
                %Fmat(idof,idof) = interfData(i).area * I3;
                indF = indF + 3;
                irowF(indF-2:indF) = idof;
                jcolF(indF-2:indF) = idof;
                valsF(indF-2:indF) = interfData(i).area;
            end
        elseif (tplas(i))
            tlim = cohes - sol_l(3*(i-1)+1)*tan(phi);
            ur = dulocTot(3*(i-1)+(1:3)) - dulocTot0(3*(i-1)+(1:3));

            % Just frictional components
            gT = ur(2:3);
            durnrm = norm(gT);
            idof = 3*(i-1)+(2:3);

            if (~(iter == 1 && tplasnew(i)) && durnrm > tol_duT)
                if (cptJ)
                    dt_dur = tlim*(durnrm^2*I2 - gT*gT')/(durnrm^3);
                    dt_dur = sparse(dt_dur);
                    Bloc = B(:,idof);

                    %Emat(idof,:) = - dt_dur*Bloc';
                    Elocal = -dt_dur*Bloc';
                    [iloc,jloc,vloc] = find(Elocal(1,:));
                    shift = length(iloc);
                    indE = indE + shift;
                    irowE(indE-shift+1:indE) = idof(1);
                    jcolE(indE-shift+1:indE) = jloc;
                    valsE(indE-shift+1:indE) = vloc;
                    [iloc,jloc,vloc] = find(Elocal(2,:));
                    shift = length(iloc);
                    indE = indE + shift;
                    irowE(indE-shift+1:indE) = idof(2);
                    jcolE(indE-shift+1:indE) = jloc;
                    valsE(indE-shift+1:indE) = vloc;

                    %Fmat(idof,3*(i-1)+1) = + interfData(i).area * tan(phi)*(gT/durnrm);
                    indF = indF + 2;
                    irowF(indF-1:indF) = idof;
                    jcolF(indF-1:indF) = 3*(i-1)+1;
                    valsF(indF-1:indF) = interfData(i).area * tan(phi)*(gT/durnrm);

                    %Fmat(idof,idof) = interfData(i).area * I2;
                    indF = indF + 2;
                    irowF(indF-1:indF) = idof;
                    jcolF(indF-1:indF) = idof;
                    valsF(indF-1:indF) = interfData(i).area;
                end
                res_l(idof(1)) = interfData(i).area * (sol_l(idof(1)) - tlim*gT(1)/durnrm);
                res_l(idof(2)) = interfData(i).area * (sol_l(idof(2)) - tlim*gT(2)/durnrm);
            else
                vaux = zeros(2,1);
                vaux(1) = sol_l(idof(1));
                vaux(2) = sol_l(idof(2));
                vauxnrm = norm(vaux);

                if (vauxnrm > 0.0)
                    if (cptJ)
                        %Fmat(idof,idof) = interfData(i).area * I2;
                        indF = indF + 2;
                        irowF(indF-1:indF) = idof;
                        jcolF(indF-1:indF) = idof;
                        valsF(indF-1:indF) = interfData(i).area;
                    end
                    res_l(idof(1)) = interfData(i).area * (sol_l(idof(1)) - tlim*vaux(1)/vauxnrm);
                    res_l(idof(2)) = interfData(i).area * (sol_l(idof(2)) - tlim*vaux(2)/vauxnrm);
                else
                    if (cptJ)
                        %Fmat(idof,idof) = interfData(i).area * I2;
                        indF = indF + 2;
                        irowF(indF-1:indF) = idof;
                        jcolF(indF-1:indF) = idof;
                        valsF(indF-1:indF) = interfData(i).area;
                    end
                    res_l(idof(1)) = 0.0;
                    res_l(idof(2)) = 0.0;
                end
            end
        end
    end

    res = [res_u;res_l];

    if (cptJ)
        % Initialize matrix with stabilization
        J = [K,B;BtDir,-C];

        Emat = sparse(irowE(1:indE),jcolE(1:indE),valsE(1:indE),3*ni,3*nn);
        Fmat = sparse(irowF(1:indF),jcolF(1:indF),valsF(1:indF),3*ni,3*ni);

        % Add Jacobian term due to sliding part
        J(3*nn+1:3*nn+3*ni,1:3*nn) = J(3*nn+1:3*nn+3*ni,1:3*nn) + Emat;
        J(3*nn+1:3*nn+3*ni,3*nn+1:3*nn+3*ni) = J(3*nn+1:3*nn+3*ni,3*nn+1:3*nn+3*ni) + Fmat;

        % Set BC on res and J
        [J,res] = DirBC(ndir,dir,zeros(ndir,1),J,res);
    else
        % Set BC on res only
        res(dir(1:ndir)) = 0;
    end

end

function [nplas,tplas,tplasnew,checkActiveSet,skipActiveSet,itA,statusL] = ...
    activeSet(Bt,indU,indL,nplas,tplas,itA,statusL,areaiR,cohes,phi,tol_duNc,tol_sig,conv,sol)
    fileID = fopen("output.txt","a");
    if (fileID == -1)
        error("Cannot open the file");
    end
    alpha = 0.05;
    skipActiveSet = false;

    % Relative displacements (w/o stabilization effects)
    dulocTot = areaiR*(Bt*sol(indU));

    ni = length(nplas);

    lambda_fix = [nplas;tplas];

    itA = itA + 1;
    statusL(:,itA) = lambda_fix;

    lambda_fix0 = lambda_fix;
    tplasnew = false(ni,1);

    sol_l = sol(indL);
    for i = 1 : ni
        duloc = dulocTot(3*(i-1)+(1:3));
        if (nplas(i))
            if (duloc(1) > -tol_duNc)
                nplas(i) = true;
                tplas(i) = true;
            else
                nplas(i) = false;
                tplas(i) = false;
            end
        elseif (sol_l(3*(i-1)+1) >= tol_sig)
            nplas(i) = true;
            tplas(i) = true;
        else
            nplas(i) = false;
            tau = sqrt(sol_l(3*(i-1)+2)^2 + sol_l(3*(i-1)+3)^2);
            tlim = cohes - sol_l(3*(i-1)+1)*tan(phi);
            if (~tplas(i) && (tau >= tlim))
                tau = (1.0-alpha) * tau;
            elseif (tplas(i) && (tau <= tlim))
                tau = (1.0+alpha) * tau;
            end
            if (tau > tlim)
                % Check if this element starts sliding
                tplasnew(i) = ~tplas(i);
                tplas(i) = true;
            else
                tplas(i) = false;
            end
        end
    end

    lambda_fix = [nplas;tplas];
    checkActiveSet = all(lambda_fix == lambda_fix0);

    if (conv)
        % Maximum allowed number of oscillations
        if (sum(ismember(statusL(:,1:itA-1)', lambda_fix', 'rows')) >= 2)
            fprintf(fileID, "       *** Active set is oscillating\n *** ");
            checkActiveSet = true;
            skipActiveSet = true;
            % Restore current (converged) solution
            nplas = lambda_fix0(1:ni);
            tplas = lambda_fix0(ni+(1:ni));
        end
    end
    fclose(fileID);

end
