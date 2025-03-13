function x = ss_newton_toy(a, x0, tol, max_iter, adaptive_choice)
    if nargin < 3
        tol = 1e-6;
    end
    if nargin < 4
        max_iter = 50;
    end
    if nargin < 5
        adaptive_choice = 'A'; % Default to Choice A
    end
    
    x = x0;
    epsilon = 1e-4; % Regularization parameter for Choice B
    
    for k = 1:max_iter
        S_x = a' * x;  % Compute S(x)
        F_x = relu(S_x) * a;  % Compute F(x)
        fprintf("iter: %d", k);
        % Compute B-subdifferential-based Jacobian based on adaptive choice
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
        
        % Check convergence
        if norm(F_x) < tol
            fprintf('Converged in %d iterations.\n', k);
            return;
        end
        
        % Check for singularity
        if rank(J_F_x) == 0
            fprintf('Singular Jacobian encountered! Stopping.\n');
            return;
        end
        
        % Solve Newton step
        delta_x = -J_F_x \ F_x;
        x = x + delta_x;
    end
    
    fprintf('Max iterations reached!\n');
end

function y = relu(x)
    y = max(0, x);
end

