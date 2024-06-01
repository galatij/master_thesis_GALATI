function dump_VTK(IDprint,nodeFieldNames,nodeScalarFields,cellFieldNames,...
                  cellScalarFields,cellFieldNames2D,cellScalarFields2D,...
                  ngauss,coord,elem,Emat,eleVol,E,nu,interfData,Bt,areaiR,...
                  sol_u,sol_l,nplas,tplas)

   % Set amplification factor for deformations 
   DEFO_AMPLI_FACTOR = 1;

   nnod = size(coord,1);
   %ne = size(elem,1);
   ni = size(interfData,1);
   ndofu = 3*nnod;
   nl = 3*ni;

   nodeFieldNames{1} = 'ux';
   nodeFieldNames{2} = 'uy';
   nodeFieldNames{3} = 'uz';
   nodeScalarFields{1} = sol_u(1:3:ndofu);
   nodeScalarFields{2} = sol_u(2:3:ndofu);
   nodeScalarFields{3} = sol_u(3:3:ndofu);
   cellFieldNames2D{1} = 'sigma';
   cellFieldNames2D{2} = 'tau1';
   cellFieldNames2D{3} = 'tau2';
   cellFieldNames2D{4} = '|tau|';
   cellFieldNames2D{5} = 'duN';
   cellFieldNames2D{6} = 'duT1';
   cellFieldNames2D{7} = 'duT2';
   cellFieldNames2D{8} = '|duT|';
   cellFieldNames2D{9} = 'nplas';
   cellFieldNames2D{10} = 'tplas';
   cellScalarFields2D{1} = sol_l(1:3:nl);
   cellScalarFields2D{2} = sol_l(2:3:nl);
   cellScalarFields2D{3} = sol_l(3:3:nl);
   cellScalarFields2D{4} = sqrt(sol_l(2:3:nl).^2+sol_l(3:3:nl).^2);

   % Compute relative displacements (global and normal)
   dulocTot = areaiR*(Bt*sol_u);% - C*sol_l);
   ur = reshape(dulocTot,3,ni)';
   cellScalarFields2D{5} = ur(:,1);
   cellScalarFields2D{6} = ur(:,2);
   cellScalarFields2D{7} = ur(:,3);
   cellScalarFields2D{8} = sqrt(ur(:,2).^2+ur(:,3).^2);
   cellScalarFields2D{9} = nplas;
   cellScalarFields2D{10} = tplas;

   % Compute stress on the elements
   stress = cpt_stress(ngauss,coord,elem,E,nu,sol_u);
   stress_fault = zeros(ni,3);
   for i = 1:ni
      top = interfData(i).etop;
      bot = interfData(i).ebottom;
      Ri = interfData(i).R_l2g;
      nrm_fault = interfData(i).normal;
      stress_loc = (stress(top,:)*eleVol(top) + stress(bot,:)*eleVol(bot))/(eleVol(top) + eleVol(bot));
      stress_tens(1,1) = stress_loc(1);
      stress_tens(1,2) = stress_loc(4);
      stress_tens(1,3) = stress_loc(5);
      stress_tens(2,1) = stress_loc(4);
      stress_tens(2,2) = stress_loc(2);
      stress_tens(2,3) = stress_loc(6);
      stress_tens(3,1) = stress_loc(5);
      stress_tens(3,2) = stress_loc(6);
      stress_tens(3,3) = stress_loc(3);
      stress_tens = Ri'*stress_tens*Ri;
      stress_fault(i,:) = (stress_tens*nrm_fault)';
   end

   cellFieldNames{1} = 'Material_ID';
   cellFieldNames{2} = 'Young_modulus';
   cellFieldNames{3} = 'SX';
   cellFieldNames{4} = 'SY';
   cellFieldNames{5} = 'SZ';
   cellFieldNames{6} = 'TXY';
   cellFieldNames{7} = 'TXZ';
   cellFieldNames{8} = 'TYZ';
   cellScalarFields{1} = Emat;
   cellScalarFields{2} = E;
   cellScalarFields{3} = stress(:,1);
   cellScalarFields{4} = stress(:,2);
   cellScalarFields{5} = stress(:,3);
   cellScalarFields{6} = stress(:,6);
   cellScalarFields{7} = stress(:,5);
   cellScalarFields{8} = stress(:,4);
   IDprint = IDprint + 1;
   fileName = sprintf('result_%3.3i.vtk', IDprint);
   coordDef = coord + DEFO_AMPLI_FACTOR*reshape(sol_u,[3,nnod])';
   write_vtk(fileName, IDprint, coordDef, elem, vertcat(interfData.top), nodeScalarFields, ...
             nodeFieldNames, cellScalarFields, cellFieldNames, ...
             cellScalarFields2D, cellFieldNames2D);

end

%-----------------------------------------------------------------------------------------

function write_vtk(fileName, time, coord3D, topol3D, topol2D, nodeScalarFields3D, ...
                   nodeFieldNames3D, cellScalarFields3D, cellFieldNames3D, ...
                   cellScalarFields2D, cellFieldNames2D)

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
    %topol2D = topol2D(:,1:4);
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

    cellScalarFields = cell(length(cellScalarFields3D)+length(cellScalarFields2D),1);
    cellFieldNames = cell(length(cellScalarFields3D)+length(cellScalarFields2D),1);
    for i = 1 : length(cellScalarFields3D)
        cellFieldNames{i} = cellFieldNames3D{i};
        if (strcmpi(cellFieldNames3D{i}, 'Material_ID'))
            cellScalarFields{i} = [cellScalarFields3D{i};ones(nFaces,1)+max(cellScalarFields3D{i})];
        else
            cellScalarFields{i} = [cellScalarFields3D{i};zeros(nFaces,1)+mean(cellScalarFields3D{i})];
        end
    end
    for i = 1 : length(cellScalarFields2D)
        cellFieldNames{i+length(cellScalarFields3D)} = cellFieldNames2D{i};
        cellScalarFields{i+length(cellScalarFields3D)} = ...
            [zeros(nCells,1)+mean(cellScalarFields2D{i});cellScalarFields2D{i}];
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
