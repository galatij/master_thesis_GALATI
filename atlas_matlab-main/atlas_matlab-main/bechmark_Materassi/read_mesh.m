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

    interfS = sort(interf,2);
    ID = zeros(ne,1);
    ID8 = 1:8;
    k = 0;
    irow = zeros(2*ni,1);
    jcol = zeros(2*ni,1);
    coef = zeros(2*ni,1);
    for i = 1 : ne
        if (matID(i) == 1)
            % Bottom
            old = topol0(i,:);
            [l1,ind] = ismember(old,list_in);
            new = old;
            new(ID8(l1)) = nn0+ind(l1);
            topol(i,:) = new;
            [l1,ind] = ismember(old,list);
            if (sum(l1) == 4)
                [~,id] = ismember(sort(topol0(i,ID8(l1))),interfS,'row');
                k = k + 1;
                irow(k) = id;
                jcol(k) = i;
                coef(k) = 2;
            end
        else
            % Top (original numbering)
            l1 = ismember(topol0(i,:),list);
            if (sum(l1) == 4)
                [~,id] = ismember(sort(topol0(i,ID8(l1))),interfS,'row');
                k = k + 1;
                irow(k) = id;
                jcol(k) = i;
                coef(k) = 1;
            end
        end
    end
    interf2e = sparse(irow,jcol,coef);

    interfB = zeros(ni,4);
    for i = 1 : ni
        [~,ind] = ismember(interf(i,:),list_in);
        interfB(i,:) = nn0+ind;
        interfB(i,ind==0) = interf(i,ind==0);
    end
    interf = [interf,interfB];

    % Dirichlet condition
    dir = zeros(3*nn+3*ni,1);
    dirval = zeros(3*nn+3*ni,1);
    I = [1:nn]';

    xmin = min(coord(:,1));
    xmax = max(coord(:,1));
    ymin = min(coord(:,2));
    ymax = max(coord(:,2));
    zmin = min(coord(:,3));
    zmax = max(coord(:,3));
    Lx = (xmax - xmin);
    Ly = (ymax - ymin);

    ndir = 0;
    % (0,[ymax-Dy,ymax])
    fac = 1000;
    ndirLoc = 0;
    while (ndirLoc == 0)
        tmp = (abs(coord(:,1)-0.0) < Lx/fac) & (coord(:,2) > ymax - Ly/20.0);
        fac = fac/2;
        ndirLoc = sum(tmp);
    end
    ndir0 = ndir;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+1;
    dirval(ndir0+1:ndir) = 0;

    % (0,[ymin,ymin+Dy])
    fac = 1000;
    ndirLoc = 0;
    while (ndirLoc == 0)
        tmp = (abs(coord(:,1)-0.0) < Lx/fac) & (coord(:,2) < ymin + Ly/20.0);
        fac = fac/2;
        ndirLoc = sum(tmp);
    end
    ndir0 = ndir;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+1;
    dirval(ndir0+1:ndir) = 0;

    % ([xmax-Dx,xmax],0)
    fac = 1000;
    ndirLoc = 0;
    while (ndirLoc == 0)
        tmp = (abs(coord(:,2)-0.0) < Ly/fac) & (coord(:,1) > xmax - Lx/20.0);
        fac = fac/2;
        ndirLoc = sum(tmp);
    end
    ndir0 = ndir;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+2;
    dirval(ndir0+1:ndir) = 0;

    % ([xmix,xmin+Dx],0)
    fac = 1000;
    ndirLoc = 0;
    while (ndirLoc == 0)
        tmp = (abs(coord(:,2)-0.0) < Ly/fac) & (coord(:,1) < xmin + Lx/20.0);
        fac = fac/2;
        ndirLoc = sum(tmp);
    end
    ndir0 = ndir;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+2;
    dirval(ndir0+1:ndir) = 0;

    % Front (fix z)
    ndir0 = ndir;
    tmp = abs(coord(:,3)-zmin) < tol;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+3;
    dirval(ndir0+1:ndir) = 0;

    % Rear (fix z)
    ndir0 = ndir;
    tmp = abs(coord(:,3)-zmax) < tol;
    ndir = ndir0 + sum(tmp);
    dir(ndir0+1:ndir) = 3*(I(tmp>0)-1)+3;
    dirval(ndir0+1:ndir) = 0;

    % just for compatibility
    dx = 0;
    dy = 0;
    dz = 0;
    coordOrig = coord;

end
