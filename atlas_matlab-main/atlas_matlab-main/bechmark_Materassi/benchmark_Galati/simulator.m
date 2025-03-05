function [itGlo, convAll] = ...
    simulator(steps, loads, bounds, nn, ni, K, B, rhs, gamma, Theta, interfData, areaiR, state0, ...
              ndir, dir, dirval, cohes, phi, noConvItMax, itmax_NR, tol_NR, ...
              maxBackStep, tol_sig, tol_duNc, tol_duT, SAVEVTK, fac, ngauss, coord, ne, ...
              topol, E, nu, volumes, matID, interf, edgeData, f2e, nrm_fault)

    indU = (1:3*nn);
    %indL = (3*nn) + (1:3*ni);
    %ntot = 3*nn + 3*ni;

    nplas0 = false(ni,1);
    tplas0 = false(ni,1);

    convAll = [];
    sol0 = state0;
    itGlo = 0;

    if (SAVEVTK)
        IDprint = 0;
        nodeFieldNames = cell(3,1);
        nodeScalarFields = cell(3,1);
        cellFieldNames = cell(8,1);
        cellScalarFields = cell(8,1);
        cellFieldNames2D = cell(10,1);
        cellScalarFields2D = cell(10,1);
    end

    for iStep = 1 : length(steps)-1
        tstart = steps(iStep);
        tend = steps(iStep+1);
        nBackStep = 0;

        dt0 = tend - tstart;
        dt = dt0;

        tcurr = tstart;
        dtcut = false;

        while (tcurr < tend)

            fprintf('t0 = %15.6e | t = %15.6e | dt = %15.6e\n', tcurr, tcurr+dt, dt);

            tcurr = tcurr + dt;

            facU = loads(indU,iStep) + (loads(indU,iStep+1)-loads(indU,iStep))*(tcurr-tstart)/dt0;
            % facL = loads(indL,iStep) + (loads(indL,iStep+1)-loads(indL,iStep))*(tcurr-tstart)/dt0;
            facDir = bounds(iStep) + (bounds(iStep+1)-bounds(iStep))*(tcurr-tstart)/dt0;            % ????????
            rhsUcurr = facU.*rhs(indU);
            % rhsLcurr = facL.*rhs(indL);
            dirvalcurr = facDir*dirval;
            rhscurr = rhsUcurr; %[rhsUcurr;rhsLcurr];

            % Solve the non linear contact mechanics problem
            % +       + +   +   +     +
            % (K - theta*B)U + N(U) + M(U) = rhs
            % +       + +   +   +     +
            % --> NOTE: output is the global solution, not just the increment
            [sol, nplas, tplas, convNRvec, convFlag, C, iter] = ...
                solve_NL_CM(nn, ni, K, theta, rhscurr, interfData, areaiR, sol0, state0, nplas0, tplas0, ...
                            ndir, dir, dirvalcurr, noConvItMax, itmax_NR, tol_NR, ...
                            cohes, phi, tol_sig, tol_duNc, tol_duT, itGlo, ngauss, coord, topol, ...
                            volumes, edgeData, f2e, E, nu);

            % backstep algo
            if (~convFlag)
                update = false;
                tcurr = tcurr - dt;
                dt = dt * 0.5;
                dtcut = true;
                nBackStep = nBackStep + 1;
            else
                update = true;
                if (dtcut)
                    dt = dt * 2.0;
                end
                dtcut = false;
            end
            
            % advance in time
            if (tcurr + dt > tend)
                dt = tend - tcurr;
            end

            if (update)     % update = true only if solve_NL_CM converged
                nplas0 = nplas;     % nplas <--> open
                tplas0 = tplas;
                itGlo = itGlo + iter;

                %sol = [sol_u;sol_l];
                sol0 = sol;

                convAll = [convAll;convNRvec];
            end

            if (nBackStep > maxBackStep)
                fprintf('Too many back steps!\n');
                return
            end

        end

        if (SAVEVTK)
            % TODO: post-process sigma
            nodeFieldNames{1} = 'ux';
            nodeFieldNames{2} = 'uy';
            nodeFieldNames{3} = 'uz';
            nodeScalarFields{1} = sol_u(1:3:3*nn);
            nodeScalarFields{2} = sol_u(2:3:3*nn);
            nodeScalarFields{3} = sol_u(3:3:3*nn);
            cellFieldNames2D{1} = 'sigma';
            cellFieldNames2D{2} = 'tau1';
            cellFieldNames2D{3} = 'tau2';
            cellFieldNames2D{4} = '|tau|';
            cellFieldNames2D{5} = 'duN';
            cellFieldNames2D{6} = 'duT1';
            cellFieldNames2D{7} = 'duT2';
            cellFieldNames2D{8} = '|duT|';
            cellFieldNames2D{9} = 'nplas';
            cellFieldNames2D{10} = 'tplas';
            cellScalarFields2D{1} = sol_l(1:3:3*ni);
            cellScalarFields2D{2} = sol_l(2:3:3*ni);
            cellScalarFields2D{3} = sol_l(3:3:3*ni);
            cellScalarFields2D{4} = sqrt(sol_l(2:3:3*ni).^2+sol_l(3:3:3*ni).^2);
            % Compute relative displacements (global and normal)
            dulocTot = areaiR*(Bt*sol(indU) - C*sol(indL));
            ur = reshape(dulocTot,3,ni)';
            cellScalarFields2D{5} = ur(:,1);
            cellScalarFields2D{6} = ur(:,2);
            cellScalarFields2D{7} = ur(:,3);
            cellScalarFields2D{8} = sqrt(ur(:,2).^2+ur(:,3).^2);
            cellScalarFields2D{9} = nplas;
            cellScalarFields2D{10} = tplas;

            stress_fault = zeros(ni,3);
            stress = cpt_stress(ngauss, coord, ne, topol, E, nu, sol_u);
            for i = 1 : ni
                top = interfData(i).etop;
                bot = interfData(i).ebottom;
                Ri = interfData(i).R;
                stress_loc = ((stress(top,:)*volumes(top) + stress(bot,:)*volumes(bot)) / (volumes(top) + volumes(bot)));
                stress_tens(1,1) = stress_loc(1);
                stress_tens(1,2) = stress_loc(4);
                stress_tens(1,3) = stress_loc(5);
                stress_tens(2,1) = stress_loc(4);
                stress_tens(2,2) = stress_loc(2);
                stress_tens(2,3) = stress_loc(6);
                stress_tens(3,1) = stress_loc(5);
                stress_tens(3,2) = stress_loc(6);
                stress_tens(3,3) = stress_loc(3);
                stress_tens = Ri'*stress_tens*Ri;
                stress_fault(i,:) = (stress_tens*nrm_fault)';
            end

            cellFieldNames{1} = 'Material_ID';
            cellFieldNames{2} = 'Young_modulus';
            cellFieldNames{3} = 'SX';
            cellFieldNames{4} = 'SY';
            cellFieldNames{5} = 'SZ';
            cellFieldNames{6} = 'TXY';
            cellFieldNames{7} = 'TXZ';
            cellFieldNames{8} = 'TYZ';
            cellScalarFields{1} = matID;
            cellScalarFields{2} = E;
            cellScalarFields{3} = stress(:,1);
            cellScalarFields{4} = stress(:,2);
            cellScalarFields{5} = stress(:,3);
            cellScalarFields{6} = stress(:,6);
            cellScalarFields{7} = stress(:,5);
            cellScalarFields{8} = stress(:,4);
            IDprint = IDprint + 1;
            fileName = sprintf('result_%3.3i.vtk', IDprint);
            coordDef = coord + fac*reshape(sol_u,[3,nn])';
            write_vtk(fileName, iStep, coordDef, topol, interf, nodeScalarFields, ...
                      nodeFieldNames, cellScalarFields, cellFieldNames, ...
                      cellScalarFields2D, cellFieldNames2D);
        end
    end

end
