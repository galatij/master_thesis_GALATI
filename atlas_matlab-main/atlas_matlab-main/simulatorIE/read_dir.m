function ind_dir = read_dir(fileIN)

ifile = fopen(fileIN);
C = textscan(fgetl(ifile),'%s');
nx = str2num(C{1}{1});
ny = str2num(C{1}{2});
nz = str2num(C{1}{3});
fclose(ifile);
tmp = dlmread(fileIN,'',1,0);
ind_dir(1:nx) = 3*(tmp(1:nx)-1)+1;
ind_dir(nx+1:nx+ny) = 3*(tmp(nx+1:nx+ny)-1)+2;
ind_dir(nx+ny+1:nx+ny+nz) = 3*(tmp(nx+ny+1:nx+ny+nz)-1)+3;

end
