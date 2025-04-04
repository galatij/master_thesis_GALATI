%% LINEAR DISPLACEMENT BENCHMARK
clear
clc
addpath('..');

%% debug cpt_normal for the computation of sigma_n
E0 = 25000;
nu = 0.25;
n = [1;2;3];
D = cpt_elas_mat(E0, nu);
A = [0.1, 0.2, 0.3; 
     0.4, 0.5, 0.6; 
     0.7, 0.8, 0.9];
u0 = [0,0,0];

u = @(x) u0 + A*x;
grad_u = @() A;
eps_u = @()[A(1,1); A(2,2); A(3,3);
            (A(1,2) + A(2,1))/2;
            (A(1,3) + A(3,1))/2;
            (A(2,3) + A(3,2))/2];
sigma = @() D*eps_u();

S = sigma();

N = cpt_normal(n);
sigma_n = N*S;
sigma_n_cmp = [ S(1) * n(1) + S(4) * n(2) + S(5) * n(3);
      S(4) * n(1) + S(2) * n(2) + S(6) * n(3);
      S(5) * n(1) + S(6) * n(2) + S(3) * n(3) ];

% disp(sigma_n - sigma_n_cmp);

%% debug cpt_assemble_B
% Idea: set D doing nothing
%       compare against the stiffness-based integral

u = @(x) x;
eps_u = [1;1;1;0;0;0];
sigma = D*eps_u;
sigma_n = n'*N*sigma;

% compute just the integral of sigma_n(phi_i)








