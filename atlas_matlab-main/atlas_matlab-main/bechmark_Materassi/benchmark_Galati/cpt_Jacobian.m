function [res,J] = cpt_Jacobian(ngauss,coord,topol,E, nu,...
                                interfData, nodePairsData, K,B,alpha,gamma, ...
                                rhs,ndir,dir,state0,sol0,areaiR, ...
                                cohes,phi,tol_duT,tol_P,iter,sol)

    cptJ = false;
    if (nargout > 1)
        cptJ = true;
    end
    nn = size(coord,1);

    dsol0 = sol0 - state0;
    stress = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,dsol0);          % nn*6
    [stress_n, stress_t] = cpt_stress_interf(stress,nodePairsData);

    %% compute the residual at iteration k
    [C0, KKT] = cpt_KKT(ngauss, coord, topol, E, nu, ...
                          interfData, nodePairsData, gamma, alpha, ...
                          dsol0, stress_n,tol_P);

    %% do the same for Friction term ...
    % [F0, FRI] = cpt_FRI(ngauss, coord, topol, E, nu, ...
    %                      interfData, nodePairsData, gamma, alpha, ...
    %                      dsol0, stress_n,tol_P);

    res = (K - alpha*B) * dsol0 + C0 - rhs; % + F0

    J = K - alpha*B + KKT; % + FRI

    [J,res] = DirBC(ndir,dir,zeros(ndir,1),J,res);
    
     %% TODO:
%     % set Dirichlet BCs
%     % return J, res


%%
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

