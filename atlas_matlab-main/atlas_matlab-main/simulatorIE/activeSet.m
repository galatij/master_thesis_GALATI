function [nplas_new,tplas_new,newly_tplas,checkActiveSet,itA,statusL] = ...
             activeSet(Bt,itA,statusL,areaiR,cohes,phi,tol_duNc,tol_sig,max_oscil,...
                       NR_converged,sol_in,nplas_in,tplas_in)

    % Increase activeSet iteration counter
    itA = itA + 1;

    % Parameter to relax tangential tension
    alpha = 0.05;
    %skipActiveSet = false;

    % Retruieve dimensions and indices
    ni = size(nplas_in,1);
    nl = size(Bt,1);
    nu = size(Bt,2);
    indU = 1:nu;
    indL = nu+1:nu+nl;

    % Relative displacements (without stabilization effects)
    dulocTot = areaiR*(Bt*sol_in(indU));

    % Get current plasticity status and store it in statusL (status history)
    lambda_fix = [nplas_in; tplas_in];
    statusL(:,itA) = lambda_fix;

    % Init indicator for newly plasticized IE
    newly_tplas = false(ni,1);
    nplas_new = false(ni,1);
    tplas_new = false(ni,1);

    sol_l = sol_in(indL);
    for i = 1:ni
       idof = 3*(i-1)+(1:3);
       duloc = dulocTot(idof);
       if (nplas_in(i))
          %**** The element was open ****
          if (duloc(1) > -tol_duNc)
             % If there is only little penetration, keep it open
             nplas_new(i) = true;
             tplas_new(i) = true;
          else
             % If there is penetration close it and rerun Newton
             nplas_new(i) = false;
             tplas_new(i) = false;
          end
       else
          %**** The element was closed ****
          if (sol_l(idof(1)) >= tol_sig)
             % If there is large enough traction open it
             nplas_new(i) = true;
             tplas_new(i) = true;
          else
             % If there is no tractionm keep it close
             nplas_new(i) = false;
             % Check tangential stress
             tau = sqrt(sol_l(idof(2))^2 + sol_l(idof(3))^2);
             tlim = cohes(i) - sol_l(idof(1))*tan(phi(i));
             if (~tplas_in(i) && (tau >= tlim))
                % If the element was closed but tau is too large, reduce a bit tau
                % to see if it was just an oscillation
                tau = (1.0-alpha) * tau;
             elseif (tplas_in(i) && (tau <= tlim))
                % If the element was open but now tau is admissible, increase a bit tau
                % to see if it was just an oscillation
                tau = (1.0+alpha) * tau;
             end
             if (tau > tlim)
                % Mark the element if it just started sliding
                newly_tplas(i) = ~tplas_in(i);
                % Mark it as slipping
                tplas_new(i) = true;
             else
                % Mark it as stick
                tplas_new(i) = false;
             end
          end
       end
    end

    % Get new plasticity status
    lambda_fix_new = [nplas_new; tplas_new];
    checkActiveSet = all(lambda_fix_new == lambda_fix);

    % If last configuration has changed but Newton has found an equilibrated solution
    % check if this configuration has already appeared, in which case force another
    % step to avoid oscillations
    if (NR_converged && ~checkActiveSet)
        % The maximum allowed number of oscillations has been reached
        if (sum(ismember(statusL(:,1:itA-1)', lambda_fix_new', 'rows')) >= max_oscil)
           % Force convergence
           checkActiveSet = true;
        end
    end

end
