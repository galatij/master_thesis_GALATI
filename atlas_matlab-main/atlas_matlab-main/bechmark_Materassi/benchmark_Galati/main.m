clc;
close all;
clear;

% NOTE
% ur = uB - uT
% normal: from T to B

% Number of Gauss points
ngauss = 2;
% Elastic parameters
E0 = 25000.e0; % MPa
nu = 0.25;
Theta = 0;  % 1 -1
gamma = 1;  % modift gamma_h
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson parameters
% ITMAX
itmax_NR = 20 * 100;
tol_NR = [1.e-7,1.e-11];
noConvItMax = 4;
maxBackStep = 4;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVEVTK = false;
SAVEVTK = true;

% Create the mesh
nx = 8;
ny = 16;
nz = 16;
Lx = 1;
Ly = 8;
Lz = 8;
[nn,coord,coordOrig,ne,topol,ni,interf,interf2e,ndir,dir,dirval,dx,dy,dz,matID] = ...
    cpt_mesh(nx,ny,nz,Lx,Ly,Lz);

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

% Prescribed stress (Neumann BCs) 
xmin = min(coordOrig(:,1));
xmax = max(coordOrig(:,1));
zmax = max(coordOrig(:,3));
ID = 1 : nf;
ID = ID(:);
bound = repmat(struct('nf', 1, 'isbound', 1, 'ID', 1, 'values', 1), 2, 1);
bound(1).isbound = (bar_fac(:,3) == zmax & bar_fac(:,1) > 0.0); % top face, positive x
bound(1).ID = ID(bound(1).isbound);                             % IDs of faces to impose stress bc
bound(1).nf = length(bound(1).ID);
bound(1).values = zeros(bound(1).nf,3) + [0,0,-3];              % impossition of stress along z
bound(2).isbound = (bar_fac(:,1) == xmax);                      % right face
bound(2).ID = ID(bound(2).isbound);
bound(2).nf = length(bound(2).ID);
bound(2).values = zeros(bound(2).nf,3) + [-2,0,0];              % imposition of stress along x (compressive)




E = zeros(ne,1) + E0;

indU = (1:3*nn);
% indL = (1:3*ni);              %%%%%%%%%%%%%%%%%%%%%%%%%

% Assemble the system
% TODO: add contribution of \theta B at the interface elements
[K,rhs,E,volumes] = assemble_K_rhs(ngauss,nn,coord,ne,topol,E,nu,e2f,faces,faceData,bound,...
                                   interf2e, interfData, Theta, gamma);


nrm_fault = [1;0;0];
%fprintf('Mesh computed\n');

tol_sig = 1.e-2 * 1e-2;
tol_duNc = 1.e-5 * 1e-2;
tol_duT = 1.e-5 * 1e-2;
cohes = 0.0;
phi = 30.0/180.0*pi;

sol_u = zeros(3*nn,1);
% state0 = sol_u; 
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

loads = zeros(ntot,length(steps));
loads(indU(1:3:3*nn),:) = repmat(loadsScal',nn,1);
loads(indU(2:3:3*nn),:) = repmat(loadsScal',nn,1);
loads(indU(3:3:3*nn),:) = repmat(loadsScalZ',nn,1);
% loads(indL,:) = 1.0;            %%%%%%%%%%%%%%%%%%%%%%%
bounds = ones(length(steps),1);         %%% bounds is an array of ones!! different from bound, used to buid the rhs

fac = 1.e2;
[itGlo, convAll] = ...
    simulator(steps, loads, bounds, nn, ni, K, rhs, gamma, Theta, interfData, areaiR, ... % state0 ...
              ndir, dir, dirval, cohes, phi, noConvItMax, itmax_NR, tol_NR, ...
              maxBackStep, tol_sig, tol_duNc, tol_duT, SAVEVTK, fac, ngauss, coord, ne, ...
              topol, E, nu, volumes, matID, interf, edgeData, f2e, nrm_fault);
