function [res,J] = cpt_Jacobian(ngauss,coord,topol,E0,nu, ...
                                interfData, nodePairsData, K,B,alpha,gamma, ...
                                rhs,ndir,dir,state0,sol0,areaiR, ...
                                cohes,phi,tol_duT,tol_P,iter,sol)

    cptJ = false;
    if (nargout > 1)
        cptJ = true;
    end
    nn = size(coord,1);

    dsol0 = sol0 - state0;
    %% compute the residual at iteration k
    [C0, KKT] = cpt_KKT(ngauss,coord,topol,E0,nu, ...
                          interfData, nodePairsData, gamma, alpha, ...
                          dsol0, 0, tol_P);

    res = (K - alpha*B) * dsol0 + C0 - rhs; 

    if (cptJ)
        % Non-linear KKT term
        indC = 0;
        irowC = zeros(144*ni,1);
        jcolC = zeros(144*ni,1);
        valsC = zeros(144*ni,1);

        % Non-linear 
        indF = 0;
        irowF = zeros(144*ni,1);
        jcolF = zeros(144*ni,1);
        valsF = zeros(144*ni,1);
    end

    % Relative displacements (w/o stabilization effects)
    %%% take the indices at the interface for top/bottom
    
    dulocTot0 = dsol0(interf(top)) - dsol0(interf(bottom));
    dulocTot = dsol(interf(top)) - dsol(interf(bottom));

    if S_x > 0
            fprintf("   S_x = %.4f, ||x|| = %.4f\n", S_x, norm(x));
            fprintf("   S_x > 0 --> choosing J_F_x = a * a'\n");
            J_F_x = a * a'; % Full rank Jacobian
        elseif S_x < 0
            fprintf("   S_x = %.4f, ||x|| = %.4f\n", S_x, norm(x));
            fprintf("   S_x < 0 --> choosing J_F_x = 0\n");
            J_F_x = zeros(length(a)); % Zero matrix
        else
            switch adaptive_choice
                case 'A' % Choosing the Largest Descent Direction
                    if norm(x) > 1e-4
                        fprintf("   S_x = %.4f, ||x|| = %.4f\n", S_x, norm(x));
                        fprintf("   -- S_x = 0,   ||x|| big enough --> choosing J_F_x = a * a'\n");
                        J_F_x = a * a';
                    else
                        fprintf("   S_x = %.4f, ||x|| = %.4f\n", S_x, norm(x));
                        fprintf("   -- S_x = 0,   ||x|| too small --> choosing J_F_x = 0\n");
                        J_F_x = zeros(length(a));
                    end
                case 'B' % Regularization Method
                    J_F_x = a * a' + epsilon * eye(length(a));
                case 'C' % Backtracking or Line Search
                    J_F_x1 = zeros(length(a));
                    J_F_x2 = a * a';
                    
                    delta_x1 = -J_F_x1 \ F_x;
                    delta_x2 = -J_F_x2 \ F_x;
                    
                    if norm(F_x + J_F_x1 * delta_x1) < norm(F_x + J_F_x2 * delta_x2)
                        J_F_x = J_F_x1;
                    else
                        J_F_x = J_F_x2;
                    end
                otherwise
                    error('Invalid adaptive choice. Use ''A'', ''B'', or ''C''.');
            end
    end

%     for i = 1 : ni
%         % if open: modify matrix and rhs s.t. lambda = lambda, i.e. no
%         % coupling
%         if (nplas(i))
%             idof = 3*(i-1)+(1:3);
%             res_l(idof(1)) = interfData(i).area * sol_l(idof(1));
%             if (cptJ)
%                 %Fmat(idof,idof) = interfData(i).area * I3;
%                 indF = indF + 3;
%                 irowF(indF-2:indF) = idof;
%                 jcolF(indF-2:indF) = idof;
%                 valsF(indF-2:indF) = interfData(i).area;            % C(i,j) will be 0 for open interfaces --> set F in the jacobian
%             end
%         % else if sliding (and not open)
%         elseif (tplas(i))
%             tlim = cohes - sol_l(3*(i-1)+1)*tan(phi);
%             % ur: increment of du between 2 subsequent time steps
%             ur = dulocTot(3*(i-1)+(1:3)) - dulocTot0(3*(i-1)+(1:3));
% 
%             % Just frictional components
%             gT = ur(2:3);
%             durnrm = norm(gT);
%             idof = 3*(i-1)+(2:3);
% 
%             if (~(iter == 1 && tplasnew(i)) && durnrm > tol_duT)    
%                 if (cptJ)
%                     dt_dur = tlim*(durnrm^2*I2 - gT*gT')/(durnrm^3);        %% missing the term: -gT*dt'/norm(gT)^2 ???
%                     dt_dur = sparse(dt_dur);
%                     Bloc = B(:,idof);
% 
%                     %Emat(idof,:) = - dt_dur*Bloc';
%                     
%                     %Fmat(idof,3*(i-1)+1) = + interfData(i).area * tan(phi)*(gT/durnrm);
%              
%                     %Fmat(idof,idof) = interfData(i).area * I2;
% 
%                 end
%                 res_l(idof(1)) = interfData(i).area * (sol_l(idof(1)) - tlim*gT(1)/durnrm);
%                 res_l(idof(2)) = interfData(i).area * (sol_l(idof(2)) - tlim*gT(2)/durnrm);
%             else
%                 vaux = zeros(2,1);
%                 vaux(1) = sol_l(idof(1));
%                 vaux(2) = sol_l(idof(2));
%                 vauxnrm = norm(vaux);
% 
%                 if (vauxnrm > 0.0)
%                     if (cptJ)
%                         %Fmat(idof,idof) = interfData(i).area * I2;
%                         indF = indF + 2;
%                         irowF(indF-1:indF) = idof;
%                         jcolF(indF-1:indF) = idof;
%                         valsF(indF-1:indF) = interfData(i).area;
%                     end
%                     res_l(idof(1)) = interfData(i).area * (sol_l(idof(1)) - tlim*vaux(1)/vauxnrm);
%                     res_l(idof(2)) = interfData(i).area * (sol_l(idof(2)) - tlim*vaux(2)/vauxnrm);
%                 else
%                     if (cptJ)
%                         %Fmat(idof,idof) = interfData(i).area * I2;
% 
%                     end
%                     res_l(idof(1)) = 0.0;
%                     res_l(idof(2)) = 0.0;
%                 end
%             end
%         end
%     end

    if (cptJ)
        % Initialize matrix with stabilization
%         J = [K,B;BtDir,-C];
% 
%         Emat = sparse(irowE(1:indE),jcolE(1:indE),valsE(1:indE),3*ni,3*nn);
%         Fmat = sparse(irowF(1:indF),jcolF(1:indF),valsF(1:indF),3*ni,3*ni);
% 
%         % Add Jacobian term due to sliding part
%         J(3*nn+1:3*nn+3*ni,1:3*nn) = J(3*nn+1:3*nn+3*ni,1:3*nn) + Emat;
%         J(3*nn+1:3*nn+3*ni,3*nn+1:3*nn+3*ni) = J(3*nn+1:3*nn+3*ni,3*nn+1:3*nn+3*ni) + Fmat;
% 
%         % Set BC on res and J
%         [J,res] = DirBC(ndir,dir,zeros(ndir,1),J,res);
%     else
%         % Set BC on res only
%         res(dir(1:ndir)) = 0;
    end

end

