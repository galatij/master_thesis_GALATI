clc;
close all;
clear;

TEST = true;

% NOTE
% ur = uB - uT
% normal: from T to B

%% Params
% Number of Gauss points
ngauss = 2;

% Elastic parameters
E0 = 1;%25000.e0; % MPa
% E = zeros(ne,1) + E0;
nu = 0.;

lambda = E0*nu/((1+nu)*(1-2*nu));
% mu = E0/2/(1+nu)

alpha = 0;  % 1 -1
gamma = 1;  % TODO: modify gamma_h

% Newton-Raphson parameters
itmax_NR = 20 * 100;
tol_NR = [1.e-7,1.e-11];
noConvItMax = 4;
maxBackStep = 4;

% Nitche's parameters
SAVEVTK = false;
tol_sig = 1.e-2 * 1e-2;
tol_duNc = 1.e-5 * 1e-2;
tol_duT = 1.e-5 * 1e-2;
tol_P = 1.e-6;
cohes = 0.0;
phi = 30.0/180.0*pi;

%% Create the mesh
nx = 16;
ny = 16;
nz = 16;
Lx = 1;
Ly = 8;
Lz = 8;

if (TEST)
    nx = 1;
    ny = 2;
    nz = 2;
    Lx = 1;
    Ly = 4;
    Lz = 4;
end

[nn,coord,coordOrig,ne,topol,ni,interf,interf2e,ndir,dir,dirval,dx,dy,dz,matID] = ...
    cpt_mesh(nx,ny,nz,Lx,Ly,Lz);

%visualize_mesh(coord, topol, interf, "test/mesh.vtk");

bar_int = zeros(ni,3);
for i = 1 : ni
    bar_int(i,:) = mean(coord(interf(i,:),:));
end
IDkeep = true(ni,1);

[nn, coord, coordOrig, ne, topol, ni, interf, interf2e, ndir, dir, dirval] = ...
    cut_fault(nn, coord, coordOrig, ne, topol, ni, interf, interf2e, ndir, dir, dirval, IDkeep);

% New barycenters
bar_int = zeros(ni,3);
for i = 1 : ni
    bar_int(i,:) = mean(coord(interf(i,:),:));
end

bar_ele = zeros(ne,3);
for i = 1 : ne
    bar_ele(i,:) = mean(coord(topol(i,:),:));
end

[nf, faces, faceData, e2f, areaf, interfData] = ...
    faceManager(ne, topol, ni, interf, coord, interf2e, bar_ele);
[nedges, edges, edgeData, f2e, interfData] = edgeManager(ni, interf, coord, interfData);
bar_fac = zeros(nf,3);
for i = 1 : nf
    bar_fac(i,:) = mean(coordOrig(faces(i,:),:));
end
bar_edg = zeros(nedges,3);
for i = 1 : nedges
    bar_edg(i,:) = mean(coordOrig(edges(i,:),:));
end

if (length(E0) == 1)
    E = E0*ones(ne,1);
else
    E = E0;
end
[nodePairsData] = nodeManager(interf, coord, interfData, E, nu);

%% Testing
if (TEST)
    fileID = fopen('test/node_manager_debug.txt', 'w');
    fprintf(fileID,"\ntopol: %d x %d\n",ne, 8);
    fprintf(fileID, 'Coordinates (coord matrix):\n');
    for i = 1:size(coord, 1)
        fprintf(fileID, '%d: %f %f %f\n', i, coord(i, :));
    end
    fprintf(fileID, '\nTopology (topol matrix):\n');
    for i = 1:size(topol, 1)
        fprintf(fileID, '%d: %d %d %d %d %d %d %d %d\n', i, topol(i, :));
    end    
    fprintf(fileID, '\nInterface nodes (interf matrix):\n');
    for i = 1:size(interf, 1)
        fprintf(fileID, '%d: %d %d %d %d %d %d %d %d\n', i, interf(i, :));
    end

    fprintf(fileID, '\nInterface Data (interfData struct):\n');
    for i = 1:length(interfData)
        fprintf(fileID, '%d:\n', i);
        fields = fieldnames(interfData);
        for j = 1:length(fields)
            value = interfData(i).(fields{j});
            if isnumeric(value)
                fprintf(fileID, '  %s: %s\n', fields{j}, mat2str(value));
            else
                fprintf(fileID, '  %s: %s\n', fields{j}, char(value));
            end
        end
        fprintf(fileID, '\n');
    end

    fprintf(fileID, '\nNode Pairs Data (nodePairsData struct):\n');
    for i = 1:length(nodePairsData)
        fprintf(fileID, '%d:\n', i);
        fields = fieldnames(nodePairsData);
        for j = 1:length(fields)
            if (fields{j} ~= "Dtop")
                value = nodePairsData(i).(fields{j});
                if isnumeric(value)
                    fprintf(fileID, '  %s: %s\n', fields{j}, mat2str(value));
                else
                    fprintf(fileID, '  %s: %s\n', fields{j}, char(value));
                end
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

%% Prescribed stress (Neumann BCs) 
xmin = min(coordOrig(:,1));
xmax = max(coordOrig(:,1));
zmax = max(coordOrig(:,3));
ID = 1 : nf;
ID = ID(:);
bound = repmat(struct('nf', 1, 'isbound', 1, 'ID', 1, 'values', 1), 2, 1);
bound(1).isbound = (bar_fac(:,3) == zmax & bar_fac(:,1) > 0.0); % top face, positive x
bound(1).ID = ID(bound(1).isbound);                             % IDs of faces to impose stress bc
bound(1).nf = length(bound(1).ID);
bound(1).values = zeros(bound(1).nf,3) + [0,0,-3];              % imposition of stress along z
bound(2).isbound = (bar_fac(:,1) == xmax);                      % right face
bound(2).ID = ID(bound(2).isbound);
bound(2).nf = length(bound(2).ID);
bound(2).values = zeros(bound(2).nf,3) + [-2,0,0];              % imposition of stress along x (compressive)

%% Assemble the system
% TODO: multply by rotation matrix ???
[K,rhs,E,volumes] = assemble_K_rhs(ngauss,nn,coord,ne,topol,E,nu,e2f,faces,faceData,bound);
[B] = assemble_B(ngauss,coord,topol,E,nu, ...
                 interf, interfData, gamma);

if (TEST)
%     % spy(K);
%     v3 = [1;2;3];
%     %spy(B);
%     interf_nod = [nodePairsData.ntop, nodePairsData.nbottom];
%     interf_nod = sort(interf_nod);
%     interf_dof = 3*(interf_nod-1)+v3;
%     interf_dof = interf_dof(:)';
%     B_full = full(B);
%     eigB = eig(B_full);
%     % disp(min(eigB));
%     sumB = sum(B, 2);
%     % disp(sumB);
% 
%     % benchmark_linear: u = @(x)x --> eps(u) = [1;1;1;0;0;0];
%     u = zeros(3*nn,1);
%     u(1:3:3*nn-2) = coord(:,1);
%     u(2:3:3*nn-1) = coord(:,2);
%     u(3:3:3*nn) = coord(:,3);
%     result1 = B*u;
%     % result2 = B_cmp*ones(3*nn,1);
%     D = cpt_elas_mat(E0,nu);
%     % disp(D)
%     % cpt_stress and cpt_stress_interf
%     [stress, eps] = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,u);
%     [s_n, s_t] = cpt_stress_interf(stress,nodePairsData);
% 
%     % benchmark periodic: u = @(x,y,z)[sin(x), sin(y), sin(z)]
%     u = zeros(3*nn,1);
%     u(1:3:3*nn-2) = sin(coord(:,2));
%     u(2:3:3*nn-1) = sin(coord(:,2));
%     u(3:3:3*nn) = sin(coord(:,3));
%     [stress, eps] = cpt_stress(ngauss,coord,topol,interfData,nodePairsData,E,nu,u);
%     [s_n, s_t] = cpt_stress_interf(stress,nodePairsData);
%     figure(1)
%     spy(K)
%     figure(2)
%     spy(B)
% 
% 
    % global assembly
    u = zeros(3*nn,1);
    u(1:3:3*nn-2) = coord(:,1);
    r = B*u;
    %disp(norm(r));
    assert(norm(B - B', 'fro') < 1e-12)     %% pass
end

nrm_fault = [1;0;0];
%fprintf('Mesh computed\n');


sol_u = zeros(3*nn,1);
state0 = sol_u; 
areai = zeros(ni,1);
for i = 1 : ni
    areai(i) = interfData(i).area;
end
areaiR = spdiags(1.0./reshape([areai';areai';areai'],3*ni,1),0,3*ni,3*ni);

steps = [0,1,2,3,4,5,6,7,8,9];
steps = steps(:);
loadsScal = zeros(length(steps),1);
loadsScal(:) = 1;
loadsScalZ = zeros(length(steps),1);
loadsScalZ(:) = [0,1,2,3,4,5,4,3,2,1];

loads = zeros(3*nn,length(steps));
loads(1:3:3*nn,:) = repmat(loadsScal',nn,1);
loads(2:3:3*nn,:) = repmat(loadsScal',nn,1);
loads(3:3:3*nn,:) = repmat(loadsScalZ',nn,1);
% loads(indL,:) = 1.0;            %%%%%%%%%%%%%%%%%%%%%%%
bounds = ones(length(steps),1);         %%% bounds is an array of ones!! different from bound, used to buid the rhs


fac = 1.e2;
[itGlo, convAll] = ...
    simulator(steps, loads, bounds, nn, ni, K, B, rhs, gamma, alpha, interfData, areaiR, state0, ...
              ndir, dir, dirval, cohes, phi, noConvItMax, itmax_NR, tol_NR, ...
              maxBackStep, tol_sig, tol_duNc, tol_duT, tol_P, SAVEVTK, fac, ngauss, coord, ne, ...
              topol, E, nu, volumes, matID, interf, nodePairsData, edgeData, f2e, nrm_fault);
