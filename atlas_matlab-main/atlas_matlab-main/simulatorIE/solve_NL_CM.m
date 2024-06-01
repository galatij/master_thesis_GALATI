function [sol_new,nplas_new,tplas_new,convNRvec,convFlag,ittot] = ...
              solve_NL_CM(control,K,Bt,areaiR,rhs_neu,ind_dir,interfData,cohes,phi,...
                          sol_init,sol_old,nplas_old,tplas_old)
%-----------------------------------------------------------------------------------------
%
% Function to solve the non-linear contact mechanics problem
% sol_init: initial solution at the beginning of the simulation
% sol_old: approximation to the solution at the beginning of this time step
% sol_curr: current approximation to the solution 
% sol_new: new approximation to the solution 
% NOTE: sol_curr and sol_old store the same quantity, a different name is used only for
%       clarity
%
%-----------------------------------------------------------------------------------------

    % Extract control parameters
    %stabKind = control.stabKind;
    maxarm = control.maxarm;
    itmax_NR = control.itmax_NR;
    tol_NR = control.tol_NR;
    actSetItMax = control.actSetItMax;
    noConvItMax = control.noConvItMax;
    tol_sig = control.tol_sig;
    tol_duNc = control.tol_duNc;
    tol_duT = control.tol_duT;
    max_oscil = control.max_oscil;

    ni = numel(interfData);

    % Store B for fast assembly
    B = Bt';

    % Compute stabilization matrix
    if control.stabKind == 3
       % Global stabilization
       Corig = cpt_global_stab(K, Bt);
       Corig = Corig - diag(sum(Corig,1));
       % Remove rows and columns corresponding to opening/sliding
       C = setStabilizationMat(Corig,nplas_old,tplas_old);
    else
       error('Analytical and Algebraic Stabilization not implemented yet');
    end

    % Init current solution
    sol_curr = sol_old;
    nplas_curr = nplas_old;
    tplas_curr = tplas_old;
    convNR = false;
    exitFlag = false;
    checkActiveSet = false;
    ittot = 0;
    itNoConv = 0;
    convNRvec = zeros(actSetItMax*itmax_NR,1);

    % Init elasTrialStep_flag to true. This is necessary to perform only one elastic
    % trialStep for this time step. Elastic means that all IE are considered closed in
    % the tangential slip. As soon as an elastic trial step is performed
    % elasTrialStep_flag is set to false in order to prevent other elastic trial steps.
    ElasticTrialStep_FLAG = true;

    % If there are no slipping elements at the beginning of this time step, this one
    % already counts as an elastic trial step, and the flag must be deactivates.
    if (sum(tplas_curr(~nplas_curr)) == 0)
       ElasticTrialStep_FLAG = false;
    end

    % Init Active set variables
    itA = 0;
    statusL = zeros(2*ni,actSetItMax);

    % This array is needed to mark tangential states change between two consecutive
    % active set iterations. This is to prevent the application of maximum plastic
    % dissipation principle to elements with null slip in the last step. It is set
    % to false at the beginning of the step
    newly_tplas = false(ni,1);

    % Acitve set loop
    while (~exitFlag)

        % Define residual computation function
        fun_res = @(sol_new,iter_NR) cpt_residual(K,B,Bt,C,rhs_neu,ind_dir,...
                                     sol_init,sol_old,areaiR,nplas_curr,tplas_curr,...
                                     newly_tplas,cohes,phi,tol_duT,interfData,...
                                     iter_NR,sol_new);

        % Solve with Newton Raphson (NOTE: nplas and tplas are freezed)
        [sol_new,iter,rnorm,convNR,resvec] = newton_solver(fun_res,sol_curr,itmax_NR,...
                                                           tol_NR,maxarm);
%@@@@@@@@@@@@@@@@
%STOPPO
%@@@@@@@@@@@

        % Update global count of NR iterations
        ittot = ittot + iter;
        convNRvec(ittot-iter+1:ittot) = resvec;

        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
        % Print total iteration count for this step
        fprintf('%6i %6i %15.6e (%5i %5i)\n', iter, ittot, rnorm, sum(nplas_curr), sum(tplas_curr));
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

        % Check plasticity status
        [nplas_new,tplas_new,newly_tplas,checkActiveSet,itA,statusL] = ...
             activeSet(Bt,itA,statusL,areaiR,cohes,phi,tol_duNc,...
                       tol_sig,max_oscil,convNR,sol_new,nplas_curr,tplas_curr);
           C = setStabilizationMat(Corig,nplas_new,tplas_new);

        if (convNR)
           % Newton has converged
           sol_curr = sol_new;
           nplas_curr = nplas_new;
           tplas_curr = tplas_new;
           % Exit only if the active set check is positive
           exitFlag = checkActiveSet;
        elseif (ElasticTrialStep_FLAG)
           % Newton did not converge in the present status ==> Try an elastic trial step
           % (close slipping IE)
           exitFlag = false;
           tplas_curr = nplas_curr;
           % Deactivate trial step (it can be tried just once)
           ElasticTrialStep_FLAG = false;
           C = setStabilizationMat(Corig,nplas_curr,tplas_curr);
        else
           % Nothing else to do, reduce time step
           exitFlag = true;
        end

        % Increase the count of non-converged Newton
        if (~convNR)
            itNoConv = itNoConv + 1;
        end

        % Check if other exit conditions have been met
        exitFlag = exitFlag || (itA > actSetItMax) || (itNoConv > noConvItMax);

    end

    % The whole procedure is converged only if both Newton and activeSet converge
    convFlag = convNR && checkActiveSet && (itA <= actSetItMax);

    % Store Newton-Raphson convergence
    convNRvec = convNRvec(1:ittot);

end
