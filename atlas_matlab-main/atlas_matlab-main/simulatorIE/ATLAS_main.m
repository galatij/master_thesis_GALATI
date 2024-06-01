clc;
close all;
clear;

% Filename
fileName = '../mesh/atlas.msh';

% Times for the simulation (must be a vertical array)
sim_times = [0 1 2 3 4]';
sim_times = [0 1]';

% Number of Gauss points
ngauss = 2;

% Elastic parameters
E0 = 25000.e0; % MPa
nu0 = 0.25;

% Fault parameters
cohes0 = 0.0;
phi0_deg = 30.0;
lam0 = [1;0;0];

% Control parameters
control.stabKind = 3;   % 1 --> analytical, 2 --> algbraic, 3 --> global
control.maxarm = 3;
control.itmax_NR = 20 * 100;
control.tol_NR = [1.e-7,1.e-11];
control.actSetItMax = 20;
control.noConvItMax = 4;
control.maxBackStep = 4;
control.tol_sig = 1.e-2 * 1e-2;
control.tol_duNc = 1.e-5 * 1e-2;
control.tol_duT = 1.e-5 * 1e-2;
control.max_oscil = 2;
control.SAVEVTK = true;

%-----------------------------------------------------------------------------------------

% Read the name of main input files
[file_coord,file_elem,file_interf,file_dir,file_neu,file_force] = read_names();

% Read coordinates
coord = read_coord(file_coord);

% Read finite elements
[elem,Emat] = read_elem(file_elem);

% Read interface elements
[interf,IEtype,IDfag,IEmat] = read_interf(file_interf);

% Read dirichlet nodes
[ind_dir] = read_dir(file_dir);

% Read Neumman conditions
neum_force = read_force(file_neu,file_force);

% Compute baricenter of IE
ni = size(interf,1);
bar_int = zeros(ni,3);
for i = 1 : ni
    bar_int(i,:) = mean(coord(interf(i,:),:));
end

% Compute baricenter of elements
ne = size(elem,1);
bar_ele = zeros(ne,3);
for i = 1 : ne
    bar_ele(i,:) = mean(coord(elem(i,:),:));
end

% Create mapping from IE to elements
i2e = map_IE_elem(interf,elem);

% Compute face data
[faces,faceData,e2f,areaf,interfData] = faceManager(elem,interf,coord,i2e,bar_ele);

% Compute edge data
[edges,edgeData,f2edg,interfData] = edgeManager(interf,coord,interfData);

% Create array for Young modulus and Poisson ratio (omogeneous material)
E = zeros(ne,1) + E0;
nu = zeros(ne,1) + nu0;

% Create array for fault properties and initial traction (omogeneous fault)
cohes = zeros(ni,1) + cohes0;
phi = zeros(ni,1) + (pi*phi0_deg/180);
lam0 = kron(ones(ni,1),lam0);

% Number of unknowns in the solid block
nn = 3*size(coord,1);

% Assemble solid block K
[K,eleVol] = assemble_K(ngauss,coord,elem,E,nu);

% Assemble interface block Bt
Bt = assemble_Bt_IE(ngauss,coord,interf,interfData);

% Store inverse of the IE area in a diagonal matrix
areai = zeros(ni,1);
for i = 1 : ni
    areai(i) = interfData(i).area;
end
areaiR = 1.0./reshape([areai';areai';areai'],3*ni,1);
areaiR = diag(sparse(areaiR));

% Perform transient non-linear simulation
[itGlo, convAll] = simulator(sim_times,control,K,Bt,areaiR,neum_force,ind_dir,ngauss,...
                             coord,elem,Emat,eleVol,E,nu,interfData,cohes,phi,lam0);

return
