function neum_force = read_force(fileNEU,fileFORCE)

ifile = fopen(fileNEU);
C = textscan(fgetl(ifile),'%s');
nx = str2num(C{1}{1});
ny = str2num(C{1}{2});
nz = str2num(C{1}{3});
fclose(ifile);
tmp = dlmread(fileNEU,'',1,0);
ind_neu(1:nx) = 3*(tmp(1:nx)-1)+1;
ind_neu(nx+1:nx+ny) = 3*(tmp(nx+1:nx+ny)-1)+2;
ind_neu(nx+ny+1:nx+ny+nz) = 3*(tmp(nx+ny+1:nx+ny+nz)-1)+3;

% Sort indices
[ind_neu,perm] = sort(ind_neu);

% Count how many time steps
nsteps = count(fileread(fileFORCE),'TIME');

% Read forces
times = zeros(nsteps,1);
forces = zeros(nx+ny+nz,nsteps);
ifile = fopen(fileFORCE,'r');
for istep = 1:nsteps
   C = textscan(fgetl(ifile),'%s');
   times(istep) = str2num(C{1}{1});
   for i = 1:nx+ny+nz
      C = textscan(fgetl(ifile),'%s');
      forces(i,istep) = str2num(C{1}{1});
   end
   % Apply permutation to forces
   forces(:,istep) = forces(perm,istep);
end
fclose(ifile);

% Store in the structure
neum_force.times = times;
neum_force.ind_neu = ind_neu;
neum_force.forces = forces;

end
