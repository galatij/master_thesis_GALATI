function [res,Jacobian] = cpt_residual(K,B,Bt,C,rhs_ext,ind_dir,sol_init,sol_old,areaiR,...
                                       nplas_curr,tplas_curr,newly_tplas,cohes,phi,...
                                       tol_duT,interfData,iter_NR,sol_curr)
%-----------------------------------------------------------------------------------------
%
% Function to compute the present residual. The residual is the "total residual" that is
% we assemble the exteranl load increment and the internal load increment from the initial
% state (beginning of the simulation where we assume internal-external force balance) to
% the present state which is unknown.
%
% sol_init: initial solution at the beginning of the simulation
% sol_old: approximation to the solution at the beginning of this time step
% sol_curr: current approximation to the solution
%
%-----------------------------------------------------------------------------------------

   % Check if Jacobian computation is needed
   cptJ = false;
   if (nargout > 1)
      cptJ = true;
   end

   % Define an auxiliary matrix
   I2 = speye(2);

   % Retrieve some dimensions
   nu = size(K,1);
   nnod = nu/3;
   ni = numel(interfData);
   nl = 3*numel(interfData);

   % Define some index set
   IDn = 3*((1:ni)-1)+1;
   indU = 1:nu;
   indL = nu+1:nu+nl;

   if (cptJ)
      % Remove constraints for open and sliding elements
      BtDir = setCouplingMat(Bt,nplas_curr,tplas_curr);
   end

   % Previous sliding (from initial state to previous time step)
   dsol_old = sol_old - sol_init;
   gT_old = Bt*dsol_old(indU);
   gT_old(IDn) = 0;

   Cl_old = C*dsol_old(indL);
   Cl_old(IDn) = 0;

   dsol_curr = sol_curr - sol_init;
   res_u = (K*dsol_curr(indU) + B*dsol_curr(indL)) - rhs_ext(indU);
   res_l = (Bt*dsol_curr(indU) - C*dsol_curr(indL) - (gT_old - Cl_old)) - rhs_ext(indL);

   if (cptJ)
      % Init indices for the Jacobian
      indE = 0;
      irowE = zeros(48*ni,1);
      jcolE = zeros(48*ni,1);
      valsE = zeros(48*ni,1);

      indF = 0;
      irowF = zeros(3*ni,1);
      jcolF = zeros(3*ni,1);
      valsF = zeros(3*ni,1);
   end

   % Relative displacements (without stabilization effects)
   dulocTot_old = areaiR*(Bt*dsol_old(indU));
   dulocTot_curr = areaiR*(Bt*dsol_curr(indU));
   sol_curr_l = sol_curr(indL);

   % Loop over IE
   for i = 1:ni
      if (nplas_curr(i))
         idof = 3*(i-1)+(1:3);
         res_l(idof(1)) = interfData(i).area * sol_curr_l(idof(1));
         res_l(idof(2)) = interfData(i).area * sol_curr_l(idof(2));
         res_l(idof(3)) = interfData(i).area * sol_curr_l(idof(3));
         if (cptJ)
            %Fmat(idof,idof) = interfData(i).area * I3;
            indF = indF + 3;
            irowF(indF-2:indF) = idof;
            jcolF(indF-2:indF) = idof;
            valsF(indF-2:indF) = interfData(i).area;
         end
      elseif (tplas_curr(i))
         % Compute Mohr-Coulomb limit tension 
         tlim = cohes(i) - sol_curr_l(3*(i-1)+1)*tan(phi(i));
         % Compute realtive displacement between previous and current step
         % (just frictional components)
         idof = 3*(i-1)+(2:3);
         gT = dulocTot_curr(idof) - dulocTot_old(idof);
         durnrm = norm(gT);

          %@@@@@@@@ IL PRIMO PREDICATO E INUTILE 
          %@@@@@@@@ QUI NON SERVE IL TPLASNEW E NON SI EFFETTUA L'ELASTIC STEP
         if (~(iter_NR == 1 && newly_tplas(i)) && durnrm > tol_duT)
            % Assembly the maximum plastic dissipation contribution to the
            % (exact) Jacobian and the total tension will follow sliding direction
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
               valsF(indF-1:indF) = interfData(i).area * tan(phi(i))*(gT/durnrm);

               %Fmat(idof,idof) = interfData(i).area * I2;
               indF = indF + 2;
               irowF(indF-1:indF) = idof;
               jcolF(indF-1:indF) = idof;
               valsF(indF-1:indF) = interfData(i).area;
            end
            % The residual is updated subdividing tension along last displacement
            % increment directions
            res_l(idof(1)) = interfData(i).area * (sol_curr_l(idof(1)) - tlim*gT(1)/durnrm);
            res_l(idof(2)) = interfData(i).area * (sol_curr_l(idof(2)) - tlim*gT(2)/durnrm);
         else
            % SIMPLIFIED TREATMENT: it is the first step or the displacement is too small,
            % hance the tangential tension is adjusted by simply scaling previous tension
            vaux = zeros(2,1);
            vaux(1) = sol_curr_l(idof(1));
            vaux(2) = sol_curr_l(idof(2));
            vauxnrm = norm(vaux);

            if (vauxnrm > 0.0)
               % Update residual with scaled tension only
               res_l(idof(1)) = interfData(i).area * (sol_curr_l(idof(1)) - tlim*vaux(1)/vauxnrm);
               res_l(idof(2)) = interfData(i).area * (sol_curr_l(idof(2)) - tlim*vaux(2)/vauxnrm);
            else
               % If the displacement is null, the residual must be null as well
               res_l(idof(1)) = 0.0;
               res_l(idof(2)) = 0.0;
            end

            if (cptJ)
               %Fmat(idof,idof) = interfData(i).area * I2;
               indF = indF + 2;
               irowF(indF-1:indF) = idof;
               jcolF(indF-1:indF) = idof;
               valsF(indF-1:indF) = interfData(i).area;
            end

         end
      end
   end

   % Compose residual
   res = [res_u;res_l];

   if (cptJ)
      % Initialize matrix with stabilization
      Jacobian = [K,B;BtDir,-C];

      % Compute contributions from maximum plastic dissipation
      Emat = sparse(irowE(1:indE),jcolE(1:indE),valsE(1:indE),3*ni,3*nnod);
      Fmat = sparse(irowF(1:indF),jcolF(1:indF),valsF(1:indF),3*ni,3*ni);

      % Add Jacobian term due to sliding part
      indU = 1:nu;
      indL = nu+1:nu+nl;
      Jacobian(indL,indU) = Jacobian(indL,indU) + Emat;
      Jacobian(indL,indL) = Jacobian(indL,indL) + Fmat;

      % Set BC on res and J
      [Jacobian,res] = DirBC(ind_dir,Jacobian,res);

   else

      % Set BC on res only
      res(ind_dir) = 0;

   end

end
