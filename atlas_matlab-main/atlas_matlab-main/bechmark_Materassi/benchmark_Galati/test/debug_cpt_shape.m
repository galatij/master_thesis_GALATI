clc; clear; close all;

% Define the nodal coordinates for a unit cube element (8 nodes)
coord = [
    -1 -1 -1;
     1 -1 -1;
     1  1 -1;
    -1  1 -1;
    -1 -1  1;
     1 -1  1;
     1  1  1;
    -1  1  1;
];
coord = [
    -1     0     0;
     0     0     0;
     0     2     0;
    -1     2     0;
    -1     0     2;
     0     0     2;
     0     2     2;
    -1     2     2];
ngauss = 2;
[nodes,~] = gausspoints(ngauss);

% Define displacement field u = [x,y,z] for each node
u = reshape(coord', [], 1);  % Flatten coord into a column vector

% Storage for strains at each node
strain_at_nodes = zeros(8,6); % Each row corresponds to strain at a node

% Compute strain at each node
for i = 1:8
%     csi = coord(i,1);
%     eta = coord(i,2);
%     theta = coord(i,3);
    csi = nodes(2);
    eta = nodes(1);
    theta = nodes(2);
    
    % Compute B at node position
    [B, detJ] = cpt_shape(coord, csi, eta, theta, 1);
    [Nloc_top] = cpt_shape_2D(coord,csi,eta,theta);
%     Nloc_matrix = repmat(Nloc_top, 8,1);
%     loc_sol = u_matrix*Nloc_top(:); 
    % Compute strain at node i
    strain_at_nodes(i, :) = (B * u)';  % Store as row vector
end

% Display computed strains at all nodes
disp('Strain at each node:');
disp(strain_at_nodes);

% Expected strain at each node: [1 1 1 0 0 0]
expected_strain = repmat([1 1 1 0 0 0], 8, 1);
disp('Expected strain at each node:');
disp(expected_strain);

% Check if computed strains match expected strain
error_norm = norm(strain_at_nodes - expected_strain, 'fro'); % Frobenius norm
if error_norm < 1e-6
    disp('✅ Strain computation at nodes is correct!');
else
    disp('❌ Strain computation is incorrect! Check cpt_shape implementation.');
end
