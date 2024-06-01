function A = cpt_finVol_stab(ngauss, coord, topol, volumes, interfData, edgeData, f2e, E, nu);

    v3 = [1;2;3];
    [colLocDof, rowLocDof] = meshgrid(v3);
    ni = length(interfData);
    nedge = size(f2e,2);
    irow = zeros(4*9*nedge,1);
    jcol = zeros(4*9*nedge,1);
    vals = zeros(4*9*nedge,1);
    indA = 1;

    for i = 1 : nedge
        list = find(f2e(:,i));
        if (length(list) == 2)
            f1 = list(1);
            f2 = list(2);

            naT = edgeData(i).n1;
            nbT = edgeData(i).n2;
            dofaT = 3*(naT-1) + v3;
            dofbT = 3*(nbT-1) + v3;

            l1a = interfData(f1).top == naT;
            l1b = interfData(f1).top == nbT;
            l2a = interfData(f2).top == naT;
            l2b = interfData(f2).top == nbT;

            naB = interfData(f1).bottom(l1a);
            nbB = interfData(f1).bottom(l1b);
            dofaB = 3*(naB-1) + v3;
            dofbB = 3*(nbB-1) + v3;

            int_area1 = cpt_area_int(ngauss, coord, topol, interfData(f1).etop, interfData(f1).top);
            A1a = int_area1(l1a);
            A1b = int_area1(l1b);

            int_area2 = cpt_area_int(ngauss, coord, topol, interfData(f2).etop, interfData(f2).top);
            A2a = int_area2(l2a);
            A2b = int_area2(l2b);

            etop1 = interfData(f1).etop;
            ebot1 = interfData(f1).ebottom;
            h1T = max(coord(topol(etop1,:),:)) - min(coord(topol(etop1,:),:));
            h1B = max(coord(topol(ebot1,:),:)) - min(coord(topol(ebot1,:),:));
            v1T = volumes(etop1);
            v1B = volumes(ebot1);

            etop2 = interfData(f2).etop;
            ebot2 = interfData(f2).ebottom;
            h2T = max(coord(topol(etop2,:),:)) - min(coord(topol(etop2,:),:));
            h2B = max(coord(topol(ebot2,:),:)) - min(coord(topol(ebot2,:),:));
            v2T = volumes(etop2);
            v2B = volumes(ebot2);

            facE1T = E(etop1)/((1+nu)*(1-2*nu));
            facE1B = E(ebot1)/((1+nu)*(1-2*nu));
            facE2T = E(etop2)/((1+nu)*(1-2*nu));
            facE2B = E(ebot2)/((1+nu)*(1-2*nu));

            R1 = interfData(f1).R;
            KlocaT = R1'*diag((facE1T*4/9*(2-3*nu)*v1T./h1T.^2).^-1)*R1;
            KlocaB = R1'*diag((facE1B*4/9*(2-3*nu)*v1B./h1B.^2).^-1)*R1;
            fac1 = A1a*A2a*(KlocaT+KlocaB);

            R2 = interfData(f2).R;
            KlocbT = R2'*diag((facE2T*4/9*(2-3*nu)*v2T./h2T.^2).^-1)*R2;
            KlocbB = R2'*diag((facE2B*4/9*(2-3*nu)*v2B./h2B.^2).^-1)*R2;
            fac2 = A1b*A2b*(KlocbT+KlocbB);

            fac = fac1 + fac2;

            % stabilization matrix assembly
            irow(indA:indA+8) = 3*(f1-1) + rowLocDof;
            jcol(indA:indA+8) = 3*(f1-1) + colLocDof;
            vals(indA:indA+8) = + fac(:);
            indA = indA + 9;

            irow(indA:indA+8) = 3*(f1-1) + rowLocDof;
            jcol(indA:indA+8) = 3*(f2-1) + colLocDof;
            vals(indA:indA+8) = - fac(:);
            indA = indA + 9;

            irow(indA:indA+8) = 3*(f2-1) + rowLocDof;
            jcol(indA:indA+8) = 3*(f1-1) + colLocDof;
            vals(indA:indA+8) = - fac(:);
            indA = indA + 9;

            irow(indA:indA+8) = 3*(f2-1) + rowLocDof;
            jcol(indA:indA+8) = 3*(f2-1) + colLocDof;
            vals(indA:indA+8) = + fac(:);
            indA = indA + 9;
        end
    end

    indA = indA - 1;
    irow = irow(1:indA);
    jcol = jcol(1:indA);
    vals = vals(1:indA);

    A = sparse(irow,jcol,vals,3*ni,3*ni);
end
