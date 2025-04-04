clear; clc;
addpath("../")

E = 1;
nu = 0.3;
D = cpt_elas_mat(E, nu);
ngauss = 2;
gamma = 10;
elem = 1;

% Simple hexahedral element
coord = [
    0, 0, 0; 
    1, 0, 0;
    1, 1, 0;
    0, 1, 0;
    0, 0, 1;
    1, 0, 1;
    1, 1, 1;
    0, 1, 1
];

topol = [1, 2, 3, 4, 5, 6, 7, 8];

% Ass: contact at x = 0
face_nodes = [1, 4, 8, 5];
normal = [-1; 0; 0];

% Compute Gloc for this face
Gloc_test = cpt_Gloc(ngauss, coord, topol, elem, face_nodes, normal, gamma, D);

disp('Computed local Gloc:')
disp(Gloc_test)

% Symmetry check
if norm(Gloc_test - Gloc_test', 'fro') < 1e-12
    disp('Gloc is symmetric');
else
    warning('Gloc is NOT symmetric');
end

% Test: force balance (sum of rows should be zero)
% idea: taking a displacement u = [c,0,0] in every node --> Gloc*u = 0,
% since no stress is applied (same for other directions) --> sum along each
% row has to be = 0!!
force_balance = sum(Gloc_test, 2);
disp('Sum of rows')
disp(force_balance)

if max(abs(force_balance)) < 1e-12
    disp('Gloc satisfies force equilibrium');
else
    warning('Gloc does NOT satisfy force equilibrium');
end
