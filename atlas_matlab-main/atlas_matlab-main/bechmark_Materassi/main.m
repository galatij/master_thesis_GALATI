clc;
close all;
clear;

% NOTE
% ur = uB - uT
% normal: from T to B
filename = "output.txt";
fileID = fopen(filename, "w");
if (fileID == -1)
    error("Cannot open the file");
end

% Number of Gauss points
ngauss = 2;
% Elastic parameters
E0 = 25000.e0; % MPa
nu = 0.25;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson parameters
% ITMAX
itmax_NR = 20 * 100;
tol_NR = [1.e-7,1.e-11];
actSetItMax = 20;
noConvItMax = 4;
maxBackStep = 4;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVEVTK = false;
SAVEVTK = true;
% Mesh deformation

%l = strsplit(fileName, '/');
%str = l{end};
%[i1,i2] = regexp(str, '\d+(\.\d+)?');
%meshID = str2num(str(i1:i2));
%[nn,coord,coordOrig,ne,topol,ni,interf,interf2e,ndir,dir,dirval,dx,dy,dz,matID] = ...
%   read_mesh(fileName);

nx = 1;
ny = 8;
nz = 8;
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

xmin = min(coordOrig(:,1));
xmax = max(coordOrig(:,1));
zmax = max(coordOrig(:,3));
ID = 1 : nf;
ID = ID(:);

neu1 = -1;
neu2 = -0.01;
bound = repmat(struct('nf', 1, 'isbound', 1, 'ID', 1, 'values', 1), 2, 1);
bound(1).isbound = (bar_fac(:,3) == zmax & bar_fac(:,1) > 0.0);
bound(1).ID = ID(bound(1).isbound);
bound(1).nf = length(bound(1).ID);
bound(1).values = zeros(bound(1).nf,3) + [0,0,neu1];
bound(2).isbound = (bar_fac(:,1) == xmax);
bound(2).ID = ID(bound(2).isbound);
bound(2).nf = length(bound(2).ID);
bound(2).values = zeros(bound(2).nf,3) + [neu2,0,0];

E = zeros(ne,1) + E0;

indU = (1:3*nn);
indL = (3*nn) + (1:3*ni);

[K,rhs,E,volumes] = assemble_K_rhs(ngauss,nn,coord,ne,topol,E,nu,e2f,faces,faceData,bound);
Bt = assemble_Bt_IE(ngauss,coord,topol,interf,interfData);
K0 = K;
[M0,ntot,rhs] = assemble_global(nn,ni,K,Bt,rhs);
lam0 = -0.00;
[lam_ini] = init_lam(ni,lam0);
nrm_fault = [1;0;0];
%fprintf('Mesh computed\n');

tol_sig = 1.e-2 * 1e-2;
tol_duNc = 1.e-5 * 1e-2;
tol_duT = 1.e-5 * 1e-2;
cohes = 0.0;
phi = 30.0/180.0*pi;

sol_u = zeros(3*nn,1);
sol_l = lam_ini;
state0 = [sol_u;sol_l];
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
loads(indL,:) = 1.0;
bounds = ones(length(steps),1);

fprintf(fileID, "===================PARAMETERS==================\n");
fprintf(fileID,"E0 = %i\n", E0);
fprintf(fileID,"nu = %.2f\n", nu);
fprintf(fileID, "[nx, ny, nz] = [%i, %i, %i]\n", nx, ny, nz);
fprintf(fileID, "[Lx, Ly, Lz] = [%i, %i, %i]\n", Lx, Ly, Lz);
fprintf(fileID, "==================NEUMANN DATA==================\n");
fprintf(fileID, "neu1 = [0,0, %.2f]\n", neu1);
fprintf(fileID, "neu2= [%.2f,0, 0]\n", neu2);
fprintf(fileID, "\n====================SIMULATION==================\n")
fclose(fileID);

fac = 1.e2;
[itGlo, convAll] = ...
    simulator(steps, loads, bounds, nn, ni, K, Bt, rhs, interfData, areaiR, state0, ...
              ndir, dir, dirval, cohes, phi, actSetItMax, noConvItMax, itmax_NR, tol_NR, ...
              maxBackStep, tol_sig, tol_duNc, tol_duT, SAVEVTK, fac, ngauss, coord, ne, ...
              topol, E, nu, volumes, matID, interf, edgeData, f2e, nrm_fault);

iter = 1:length(convAll);
figure()
plot(iter, convAll);
title("Convergence History")
xlabel("it")
ylabel("res")