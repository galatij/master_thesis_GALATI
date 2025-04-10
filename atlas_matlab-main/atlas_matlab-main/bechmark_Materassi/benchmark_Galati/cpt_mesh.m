function [nn,coord,coordOrig,ne,topol,ni,interf,interf2e,ndir,dir,dirval,dx,dy,dz,matID] = ...
    cpt_mesh(nx,ny,nz,Lx,Ly,Lz)

    nnx = 2*(nx+1);     % +1 due to node duplication at fault
    nny = ny+1;
    nnz = nz+1;
    nn = nnx*nny*nnz;
    ne = 2*nx*ny*nz;
    ne_tot = (2*nx+1)*ny*nz;

    dx = Lx / nx;
    dy = Ly / ny;
    dz = Lz / nz;

    coord = zeros(nn,3);
    kk = 0;
    z = 0;
    for k = 1 : nnz
        y = 0;
        for j = 1 : nny
            x = -nx*dx;
            for i = 1 : nnx
                kk = kk + 1;
                coord(kk,1) = x;
                coord(kk,2) = y;
                coord(kk,3) = z;
                x = x + dx;
                if (i == nx+1)
                    x = x - dx; % same coord for duplicate nodes
                end
            end
            y = y + dy;
        end
        z = z + dz;
    end

    coordOrig = coord;

    topol_all = zeros(ne_tot,8);
    id = zeros(ne_tot,1);
    nxy = nnx*nny;
    for k = 1 : nz
        for j = 1 : ny
            for i = 1 : 2*nx+1
                kk = (k-1)*(2*nx+1)*ny + (j-1)*(2*nx+1) + i;
                if (i == nx+1)
                    id(kk) = 1;
                end
                i1 = (k-1)*nxy + (j-1)*nnx + i;
                topol_all(kk,1) = i1;
                topol_all(kk,2) = i1+1;
                topol_all(kk,3) = i1+1+nnx;
                topol_all(kk,4) = i1+nnx;
                topol_all(kk,5) = i1+nxy;
                topol_all(kk,6) = i1+1+nxy;
                topol_all(kk,7) = i1+1+nnx+nxy;
                topol_all(kk,8) = i1+nnx+nxy;
            end
        end
    end

    topol = topol_all(id==0,:);
    interf = topol_all(id==1,:);

    matID = zeros(ne,1);
    for i = 1 : ne
        if (mean(coordOrig(topol(i,:),1)) < 0.0)
            matID(i) = 1;
        else
            matID(i) = 2;
        end
    end

    perm = [1,4,8,5,2,3,7,6];
    interf = interf(:,perm);
    ni = size(interf,1);

    % Fix new coordinates on fault
    for i = 1 : ni
        coord(interf(i,5:8),:) = coord(interf(i,1:4),:);
    end

    % Find elements touching the fault
    irow = zeros(2*ni,1);
    jcol = zeros(2*ni,1);
    coef = zeros(2*ni,1);
    kk = 0;
    for i = 1 : ni
        j = mod(i-1,ny)+1;
        k = floor((i-1)/ny)+1;
        kk = kk + 1;
        % Top
        irow(kk) = i;
        jcol(kk) = (k-1)*2*nx*ny + (j-1)*2*nx + nx;
        coef(kk) = 1;
        kk = kk + 1;
        % Bottom
        irow(kk) = i;
        jcol(kk) = (k-1)*2*nx*ny + (j-1)*2*nx + nx + 1;
        coef(kk) = 2;
    end
    interf2e = sparse(irow, jcol, coef);

    %% Dirichlet BCs
    dir = zeros(3*nn+3*ni+ni,1); % just preallocating, not all entries are used
    dirval = zeros(3*nn+3*ni+ni,1);
    j = 0;
    % Bottom/top: fully constrain bottom face only
    for i = 1 : nnx*nny
        k = i;
        j = j + 1;
        dir(j) = 3*(k-1)+1;
        dirval(j) = 0.0;
        j = j + 1;
        dir(j) = 3*(k-1)+2;
        dirval(j) = 0.0;
        j = j + 1;
        dir(j) = 3*(k-1)+3;
        dirval(j) = 0.0;

        %k = i + nnx*nny*(nnz-1);
        %j = j + 1;
        %dir(j) = 3*(k-1)+1;
        %dirval(j) = 0.0;
        %j = j + 1;
        %dir(j) = 3*(k-1)+2;
        %dirval(j) = 0.0;
%        j = j + 1;
%        dir(j) = 3*(k-1)+3;
%        dirval(j) = 0.0;
    end
    % Left/right: constrain x direction for left only
    for i = 1 : nny*nnz
        k = (i-1)*nnx + 1;
        j = j + 1;
        dir(j) = 3*(k-1)+1;
        dirval(j) = 0.0;
        %j = j + 1;
        %dir(j) = 3*(k-1)+2;
        %dirval(j) = 0.0;
        %j = j + 1;
        %dir(j) = 3*(k-1)+3;
        %dirval(j) = 0.0;

% THESE ROWS CONSTRAINT THE SURFACE WITH NORMAL x
%        k = (i-1)*nnx + nnx;
%        j = j + 1;
%        dir(j) = 3*(k-1)+1;
%        dirval(j) = 0.0;
        %j = j + 1;
        %dir(j) = 3*(k-1)+2;
        %dirval(j) = 0.0;
        %j = j + 1;
        %dir(j) = 3*(k-1)+3;
        %dirval(j) = 0.0;
    end
    % Front/back: constrain y direction of both front and back
    for i = 1 : nnz
        k = (i-1)*nnx*nny;
        for i1 = 1 : nnx
            k1 = k + i1;
            %j = j + 1;
            %dir(j) = 3*(k1-1)+1;
            %dirval(j) = 0.0;
            j = j + 1;
            dir(j) = 3*(k1-1)+2;
            dirval(j) = 0.0;
            %j = j + 1;
            %dir(j) = 3*(k1-1)+3;
            %dirval(j) = 0.0;
        end
        for i1 = 1 : nnx
            k1 = k + (nny-1)*nnx + i1;
            %j = j + 1;
            %dir(j) = 3*(k1-1)+1;
            %dirval(j) = 0.0;
            j = j + 1;
            dir(j) = 3*(k1-1)+2;
            dirval(j) = 0.0;
            %j = j + 1;
            %dir(j) = 3*(k1-1)+3;
            %dirval(j) = 0.0;
        end
    end
    ndir = j; % gives the relevant entries in dir and dirval (number of constrained DOFs)

end
