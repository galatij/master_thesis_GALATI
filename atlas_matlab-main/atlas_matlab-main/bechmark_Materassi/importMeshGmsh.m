function mesh = importMeshGmsh(fileName)

    fid = fopen(fileName, 'r');
    cont = ~feof(fid);
    while (cont)
        str = fgetl(fid);
        cont = isempty(strfind(str, 'Nodes'));
    end
    str = fgetl(fid);
    nn = sscanf(str, '%i', 1);
    coord = fscanf(fid, '%e\n', [4,nn])';
    mesh.coord = coord(:,2:4);
    cont = ~feof(fid);
    while (cont)
        str = fgetl(fid);
        cont = isempty(strfind(str, 'Elements'));
    end
    str = fgetl(fid);
    ne = sscanf(str, '%i', 1);
    topol1D = zeros(ne,3);
    topol2D = zeros(ne,5);
    topol3D = zeros(ne,9);
    k1D = 0;
    k2D = 0;
    k3D = 0;
    for i = 1 : ne
        str = fgetl(fid);
        j = sscanf(str, '%*i %i', 1);
        if (j == 1)
            k1D = k1D + 1;
            topol1D(k1D,:) = sscanf(str, '%*i %*i %*i %i %*i %i %i', 3);
        end
        if (j == 3)
            k2D = k2D + 1;
            topol2D(k2D,:) = sscanf(str, '%*i %*i %*i %i %*i %i %i %i %i', 5);
        end
        if (j == 5)
            k3D = k3D + 1;
            topol3D(k3D,:) = sscanf(str, '%*i %*i %*i %i %*i %i %i %i %i %i %i %i %i', 9);
        end
    end
    mesh.topol1D = topol1D(1:k1D,:);
    mesh.topol2D = topol2D(1:k2D,:);
    mesh.topol3D = topol3D(1:k3D,:);
    fclose(fid);

end
