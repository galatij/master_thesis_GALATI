function write_vtk(fileName, time, coord3D, topol3D, topol2D, nodeScalarFields3D, ...
                   nodeFieldNames3D, cellScalarFields3D, cellFieldNames3D)

    % Open output file
    FID_vtk = fopen(fileName, 'wb', 'ieee-be');

    % Header
    fprintf(FID_vtk, '%s\n', '# vtk DataFile Version 3.0');
    fprintf(FID_vtk, '%s\n', 'RESULT');
    fprintf(FID_vtk, '%s\n', 'BINARY');
    fprintf(FID_vtk, '%s\n', 'DATASET UNSTRUCTURED_GRID');

    fprintf(FID_vtk, '%s\n', 'FIELD FieldData 1');
    fprintf(FID_vtk, '%s\n', 'TIME 1 1 float');
    fwrite(FID_vtk, time, 'float32');

    % Point coordinates
    nNodes = size(coord3D,1);
    fprintf(FID_vtk, 'POINTS %d float\n', nNodes);
    fwrite(FID_vtk, coord3D', 'float32');

    % Mesh topology
    nCells = size(topol3D,1);
    topol2D = topol2D(:,1:4);
    nFaces = size(topol2D,1);
    fprintf(FID_vtk, 'CELLS %d %d\n', nCells+nFaces,nCells+nCells*size(topol3D,2)+nFaces+nFaces*size(topol2D,2));
    tmp = [8*ones(1,nCells);topol3D'-1];
    fwrite(FID_vtk, uint32(tmp), 'uint32');
    tmp = [4*ones(1,nFaces);topol2D'-1];
    fwrite(FID_vtk, uint32(tmp), 'uint32');

    fprintf(FID_vtk, 'CELL_TYPES %d\n', nCells+nFaces);
    fwrite(FID_vtk, uint32(12*ones(nCells,1)), 'uint32');
    fwrite(FID_vtk, uint32(9*ones(nFaces,1)), 'uint32');

    write_attributes(FID_vtk, 'POINT_DATA', nNodes, nodeScalarFields3D, nodeFieldNames3D);

    cellScalarFields = cell(length(cellScalarFields3D),1);
    cellFieldNames = cell(length(cellScalarFields3D),1);
    for i = 1 : length(cellScalarFields3D)
        cellFieldNames{i} = cellFieldNames3D{i};
        if (strcmpi(cellFieldNames3D{i}, 'Material_ID'))
            cellScalarFields{i} = [cellScalarFields3D{i};ones(nFaces,1)+max(cellScalarFields3D{i})];
        else
            cellScalarFields{i} = [cellScalarFields3D{i};zeros(nFaces,1)+mean(cellScalarFields3D{i})];
        end
    end
    write_attributes(FID_vtk, 'CELL_DATA', nCells+nFaces, cellScalarFields, cellFieldNames);
    fclose(FID_vtk);

end

function write_attributes(FID, dataset_name, num_pts, scalarFields, fieldNames)

    fprintf(FID, '%s %d\n', dataset_name, num_pts);

    for iField = 1 : numel(fieldNames)
        fprintf(FID, 'SCALARS %s float 1\n',fieldNames{iField});
        fprintf(FID, 'LOOKUP_TABLE default\n');
        fwrite(FID, scalarFields{iField}, 'float32');
        fprintf(FID, '\n');
    end
end
