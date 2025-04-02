function visualize_mesh(coord, topol, interf, filename)

        nn = size(coord,1);
        ne = size(topol,1);
        ni = size(interf,1);

        nodeFieldNames = cell(4,1);
        nodeScalarFields = cell(3,1);
        nodeFieldNames{1} = 'x';
        nodeFieldNames{2} = 'y';
        nodeFieldNames{3} = 'z';
        nodeFieldNames{4} = 'nodeID';
        nodeScalarFields{1} = coord(:,1);
        nodeScalarFields{2} = coord(:,2);
        nodeScalarFields{3} = coord(:,3);
        nodeScalarFields{4} = [1:nn]';

        cellFieldNames = cell(1,1);
        cellScalarFields = cell(1,1);
        cellFieldNames{1} = 'elemID';
        cellScalarFields{1} = [1:ne]';

        cellFieldNames2D = cell(10,1);
        cellScalarFields2D = cell(10,1);
        cellFieldNames2D{1} = 'interfID';
        cellScalarFields2D{1} = [1:ni]';

        % coordDef = coord + fac*reshape(sol_u,[3,nn])';
        write_vtk(filename, 0, coord, topol, interf, nodeScalarFields, ...
                  nodeFieldNames, cellScalarFields, cellFieldNames, ...
                  cellScalarFields2D, cellFieldNames2D);
