% Example usage
a = [1; 2];  % Weight vector
x0 = [-1; -1];  % Initial guess
solution_A = ss_newton_toy(a, x0, 1e-6, 50, 'A');
disp('Solution using Adaptive Choice A:'), disp(solution_A);
solution_B = ss_newton_toy(a, x0, 1e-6, 50, 'B');
disp('Solution using Adaptive Choice B:'), disp(solution_B);
solution_C = ss_newton_toy(a, x0, 1e-6, 50, 'C');
disp('Solution using Adaptive Choice C:'), disp(solution_C);