function [nn,coord,coordOrig,ne,topol,ni,interf,interf2e,ndir,dir,dirval,dx,dy,dz,matID] = ...
    read_mesh(fileName)

    mesh = importMeshGmsh(fileName);

    tol = 1.e-10;

    coord0 = mesh.coord;
    interf = mesh.topol2D;
    topol0 = mesh.topol3D;
    nn0 = size(coord0,1);
    ne = size(topol0,1);

    matID = topol0(:,1);
    topol0 = topol0(:,2:9);

    bound = mesh.topol1D(:,2:3);
    bound = unique(bound(:));

    interf = interf(:,2:5);
    ni = size(interf,1);
    list = unique(interf(:));
    list_in = setdiff(list, bound);

    coord = [coord0;coord0(list_in,:)];
    nn = size(coord,1);
    topol = topol0;

    fileInterf = '../mesh/interf';
    fid = fopen(fileInterf, 'r');
    str = fgetl(fid);
    ni = sscanf(str, '%i', 1);
    interf = zeros(ni,8);
    for i = 1 : ni
        str = fgetl(fid);
        interf(i,:) = sscanf(str, '%*i %i %i %i %i %i %i %i %i', 8);
    end
    fclose(fid);

    face = [1,2,3,4;1,2,6,5;2,3,7,6;4,3,7,8;1,4,8,5;5,6,7,8];
    irow = zeros(2*ni,1);
    jcol = zeros(2*ni,1);
    coef = zeros(2*ni,1);
    intT = sort(interf(:,5:8),2);
    intB = sort(interf(:,1:4),2);
    k = 0;
    for i = 1 : ne
        for j = 1 : 6
            f = sort(topol(i,face(j,:)));
            [flag, ii] = ismember(f, intT, 'rows');
            if (flag)
                k = k + 1;
                irow(k) = ii;
                jcol(k) = i;
                coef(k) = 1;
            end
            [flag, ii] = ismember(f, intB, 'rows');
            if (flag)
                k = k + 1;
                irow(k) = ii;
                jcol(k) = i;
                coef(k) = 2;
            end
        end
    end
    interf2e = sparse(irow,jcol,coef);

    % Dirichlet condition
    fileInterf = '../mesh/dir_ind';
    fid = fopen(fileInterf, 'r');
    str = fgetl(fid);
    ndir = sscanf(str, '%i %i %i', 3);
    dirX = fscanf(fid, '%i', ndir(1));
    dirY = fscanf(fid, '%i', ndir(2));
    dirZ = fscanf(fid, '%i', ndir(3));
    fclose(fid);
    ndir = sum(ndir);
    dir = [3*(dirX-1)+1; 3*(dirY-1)+2; 3*(dirZ-1)+3];
    dirval = zeros(3*nn,1);

    % just for compatibility
    dx = 0;
    dy = 0;
    dz = 0;
    coordOrig = coord;

end
