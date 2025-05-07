function [itGlo, convAll] = ...
    simulator(steps, loads, bounds, nn, ni, K, B, rhs, gamma, alpha, interfData, areaiR, state0, ...
              ndir, dir, dirval, cohes, phi, noConvItMax, itmax_NR, tol_NR, ...
              maxBackStep, tol_sig, tol_duNc, tol_duT, tol_P, SAVEVTK, fac, ngauss, coord, ne, ...
              topol, E, nu, volumes, matID, interf, nodePairsData, nrm_fault)

    TEST = true;
    indU = (1:3*nn);

    convAll = [];
    sol0 = state0;
    itGlo = 0;

    if (SAVEVTK)
        IDprint = 0;
        nodeFieldNames = cell(15,1);
        nodeScalarFields = cell(15,1);
        cellFieldNames = cell(2,1);
        cellScalarFields = cell(2,1);
%         cellFieldNames2D = cell(10,1);
%         cellScalarFields2D = cell(10,1);
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
            facDir = bounds(iStep) + (bounds(iStep+1)-bounds(iStep))*(tcurr-tstart)/dt0;            % ????????
            rhscurr = facU.*rhs(indU);
            dirvalcurr = facDir*dirval;

            % Solve the non linear contact mechanics problem
            % 
            % (K - theta*B)U + C(U) + E(U) = rhs
            % 
            % --> NOTE: output is the global solution, not just the increment
            [sol, convNRvec, convFlag, iter] = ...
                solve_NL_CM(K, B, alpha, gamma, E, nu, rhscurr, interfData, nodePairsData, areaiR, sol0, state0, ...
                            ndir, dir, dirvalcurr, noConvItMax, itmax_NR, tol_NR, ...
                            cohes, phi, tol_sig, tol_duNc, tol_duT, tol_P, ngauss, coord, topol, ...
                            volumes);

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

            if (update)
                % nplas and tplas should now be referre to the nodes!!
                itGlo = itGlo + iter;
                sol0 = sol;
                convAll = [convAll;convNRvec];
            end

            if (nBackStep > maxBackStep)
                fprintf('Too many back steps!\n');
                return
            end

        end

        % TODO: modify post-process for output
        if (SAVEVTK)
            nodeFieldNames{1} = 'ux';
            nodeFieldNames{2} = 'uy';
            nodeFieldNames{3} = 'uz';
            nodeScalarFields{1} = sol(1:3:3*nn);
            nodeScalarFields{2} = sol(2:3:3*nn);
            nodeScalarFields{3} = sol(3:3:3*nn);

            % TODO: cpt_stress in all the nodes --> cpt traction
            stress = cpt_stress_tot(ngauss,coord,topol,E,nu,sol);
%             traction = cpt_traction_tot(. . .);
            nodeFieldNames{4} = 'Sx';
            nodeFieldNames{5} = 'Sy';
            nodeFieldNames{6} = 'Sz';
            nodeFieldNames{7} = 'Tyz';
            nodeFieldNames{8} = 'Txz';
            nodeFieldNames{9} = 'Txy';
            nodeScalarFields{4} = stress(:,1);
            nodeScalarFields{5} = stress(:,2);
            nodeScalarFields{6} = stress(:,3);
            nodeScalarFields{7} = stress(:,6);
            nodeScalarFields{8} = stress(:,5);
            nodeScalarFields{9} = stress(:,4);

            % Fault data
            nodeFieldNames{10} = 'sigma_n';
            nodeFieldNames{11} = '|tau|';
            nodeFieldNames{12} = 'duN';
            nodeFieldNames{13} = 'duT1';
            nodeFieldNames{14} = 'duT2';
            nodeFieldNames{15} = '|duT|';
            nodeScalarFields{10} = nan(nn,1);
            nodeScalarFields{11} = nan(nn,1);
            nodeScalarFields{12} = nan(nn,1);
            nodeScalarFields{13} = nan(nn,1);
            nodeScalarFields{14} = nan(nn,1);
            nodeScalarFields{15} = nan(nn,1);
            v3 = [1;2;3];
            nni = numel(nodePairsData);
            stress_interf = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,sol);
            [sigma_n, sigma_t] = cpt_stress_interf(stress_interf, nodePairsData);
            for i = 1:nni
                n = nodePairsData(i).normal;
                t1 = nodePairsData(i).t1;
                t2 = nodePairsData(i).t2;
                ntop = nodePairsData(i).ntop;
                nbot = nodePairsData(i).nbottom;
                dof_top = 3*(ntop-1)+v3;
                dof_bot = 3*(nbot-1)+v3;
                du = sol(dof_top(:)) - sol(dof_bot(:));
                nodeScalarFields{10}(ntop) = sigma_n(i);
                nodeScalarFields{11}(ntop) = norm(sigma_t(i,:));
                nodeScalarFields{12}(ntop) = n'*du;
                nodeScalarFields{13}(ntop) = t1'*du;
                nodeScalarFields{14}(ntop) = t2'*du;
                nodeScalarFields{15}(ntop) = norm(nodeScalarFields{13}(ntop)...
                                                 +nodeScalarFields{14}(ntop));
            end                

            cellFieldNames{1} = 'Material_ID';
            cellFieldNames{2} = 'Young_modulus';
            cellScalarFields{1} = matID;
            cellScalarFields{2} = E;

            IDprint = IDprint + 1;
            fileName = sprintf('result_%3.3i.vtk', IDprint);
            coordDef = coord + fac*reshape(sol,[3,nn])';
            write_vtk(fileName, iStep, coordDef, topol, interf, nodeScalarFields, ...
                      nodeFieldNames, cellScalarFields, cellFieldNames);                      
        end
    end

end
