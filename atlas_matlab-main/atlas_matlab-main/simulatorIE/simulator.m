function [itGlo, convAll] = simulator(sim_times,control,K,Bt,areaiR,neum_force,ind_dir,...
                                      ngauss,coord,elem,Emat,eleVol,E,nu,interfData,...
                                      cohes,phi,lam0)
%-----------------------------------------------------------------------------------------
%
%  The simulator executes several time steps.
%  sol_init: denotes the solution at the beginning of the simulation
%  sol_old: denotes the solution at the beginning of the current time step
%  sol_new: denotes the solution at the end of the current time step, as is computed by the
%           Newton-Raphson scheme
%
%-----------------------------------------------------------------------------------------

   % Recover main dimensions
   nnod = size(coord,1);
   %ne = size(elem,1);
   ni = size(interfData,1);
   ndofu = 3*nnod;
   nl = 3*ni;
   nn = ndofu + nl;

   % Init matrix to store convergence history
   itGlo = 0;
   convAll = [];

   if (control.SAVEVTK)
      % Init the VTK variables 
      IDprint = 0;
      nodeFieldNames = cell(3,1);
      nodeScalarFields = cell(3,1);
      cellFieldNames = cell(8,1);
      cellScalarFields = cell(8,1);
      cellFieldNames2D = cell(10,1);
      cellScalarFields2D = cell(10,1);
   end

   % Set initial state
   %time_0 = sim_times(1);
   sol_init = [zeros(ndofu,1); lam0];
   sol_old = sol_init;

   % Init open/close and stick/slip indicators
   nplas_old = false(ni,1);
   tplas_old = false(ni,1);

   % Set initial time
   time_0 = sim_times(1);

   % Compute Neumman part of rhs at time_0 (initial time)
   rhs0_neu = assemble_rhs_neu(nn,time_0,neum_force);

   % Time step loop
   T_end = time_0;
   for iStep = 2:length(sim_times)

      % Get times
      T_start = T_end;
      T_end = sim_times(iStep);
      nBackStep = 0;

      DT0 = T_end - T_start;
      DT = DT0;

      % Set current time
      T_curr = T_start;
      BackStep_flag = false;

      % Backstep loop
      while (T_curr < T_end)

         % Set the next time
         T_next = T_curr + DT;

         % Print times
         fprintf('T_curr = %15.6e | T_next = %15.6e | DT = %15.6e\n',...
                 T_curr,T_next, DT);

         % Compute increment of external load with respect to initial load
         rhs_neu = assemble_rhs_neu(nn,T_next,neum_force) - rhs0_neu;

         % Solve the non linear contact mechanics problem
         % +       + +   +   +     +
         % |  K  B | | u | = | r_f |
         % | B^T C | | l |   |  0  |
         % +       + +   +   +     +
         % --> NOTE: output is the global solution, not just the increment
         [sol_new,nplas_new,tplas_new,convNRvec,convFlag,iter] = ...
             solve_NL_CM(control,K,Bt,areaiR,rhs_neu,ind_dir,interfData,cohes,phi,...
                         sol_init,sol_old,nplas_old,tplas_old);

         % Check Newton-Raphson convergence
         if (~convFlag)
            % Fail
            update_flag = false;
            DT = DT * 0.5;
            BackStep_flag = true;
            nBackStep = nBackStep + 1;
         else
            % Success
            update_flag = true;
            if (BackStep_flag)
                DT = DT * 2.0;
            end
            BackStep_flag = false;
            T_curr = T_next;
         end

         % Control that time step do not overcome T_end
         if (T_curr + DT > T_end)
             DT = T_end - T_curr;
         end

         if (update_flag)
            % Global iteration count
            itGlo = itGlo + iter;
            % Open/close and stick/slip indicators
            nplas_old = nplas_new;
            tplas_old = tplas_new;
            % Current solution
            sol_old = sol_new;
            % Store convergence history
            convAll = [convAll;convNRvec];
         end

         if (nBackStep > control.maxBackStep)
            fprintf('Too many back steps!\n');
            return
         end

      end % End of backstep loop

      % Dump on VTK file if required
      if (control.SAVEVTK)
          sol_u = sol_new(1:ndofu);
          sol_l = sol_new(ndofu+1:end);
          dump_VTK(IDprint,nodeFieldNames,nodeScalarFields,cellFieldNames,...
                  cellScalarFields,cellFieldNames2D,cellScalarFields2D,...
                  ngauss,coord,elem,Emat,eleVol,E,nu,interfData,Bt,areaiR,...
                  sol_u,sol_l,nplas_new,tplas_new);
      end

   end % End of time step loop

end
