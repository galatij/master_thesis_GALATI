function [nn1, coord1, coordOrig1, ne, topol1, ni1, interf1, interf2e1, ndir1, dir1, dirval1] = ...
    cut_fault(nn, coord, coordOrig, ne, topol, ni, interf, interf2e, ndir0, dir, dirval, IDkeep);

    interf2e1 = interf2e(IDkeep,:);
    ni1 = sum(IDkeep);

    topol1 = topol;
    map1 = zeros(nn,1);
    map1(interf(~IDkeep,5:8)) = interf(~IDkeep,1:4);
    map2 = zeros(nn,1);

    j = 0;
    for i = 1 : nn
        if (map1(i) == 0)
            j = j + 1;
            map2(i) = j;
        end
    end
    nn1 = j;

    coord1 = coord(map2>0,:);
    coordOrig1 = coordOrig(map2>0,:);

    map2_3D = 3*(map2-1)+[1,2,3];
    map2_3D(map2_3D<0) = 0;
    map2_3D = map2_3D';
    map2_3D = map2_3D(:);
    dir_loc = zeros(3*nn,1);
    dir_loc(dir(dir>0)) = 1;
    tmp = dir_loc(map2_3D>0);
    ID = 1 : length(tmp);
    ID = ID(:);
    dir1 = ID(tmp>0);
    dirval_loc = zeros(3*nn,1);
    dirval_loc(dir(dir>0)) = dirval(dir>0);
    dirval1 = dirval_loc(tmp>0);
    ndir1 = length(dir1);
    dir1 = [dir1;zeros(3*nn1+3*ni1+ni1-ndir1,1)];
    dirval1 = [dirval1;zeros(3*nn1+3*ni1+ni1-ndir1,1)];

    map1(map1>0) = map2(map1(map1>0));
    map2 = map1 + map2;
    topol1 = map2(topol);
    interf1 = map2(interf);
    interf1 = interf1(IDkeep,:);

    map2_3D = 3*(map2-1)+[1,2,3];
    map2_3D(map2_3D<0) = 0;
    map2_3D = map2_3D';
    map2_3D = map2_3D(:);

end
