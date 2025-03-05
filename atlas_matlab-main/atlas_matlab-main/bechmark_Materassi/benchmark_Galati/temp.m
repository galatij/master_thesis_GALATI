function [B] = assemble_B(ngauss,coord,ne,topol,E0,nu,e2f,faces,faceData,bound, interf, interfData);

    [nodes,weights] = gausspoints(ngauss);
    ni = size(interf,1);
    nn = size(coord,1);
    if (length(E0) == 1)
        E = E0*ones(ne,1);
    else
        E = E0;
    end

    f2e = e2f';

    k = 1;
    Blist = zeros(ne*576,3);
    v3 = [1;2;3];
    face = [1,2,3,4;1,2,6,5;2,3,7,6;4,3,7,8;1,4,8,5;5,6,7,8];
    volumes = zeros(ne,1);

    for i = 1 : ne
        % Find elements corresponding to the interface 
        loc_nod = topol(i,:);
        loc_faces = e2f(loc_nod);

        % If the face is at the interface


        % Local stiffness matrix
        D = cpt_elas_mat(E(i),nu);


        loc_coo = coord(loc_nod,:);
        Bloc = zeros(24,24);
        vol = 0.0;
        for i1 = 1 : ngauss
            csi = nodes(i1);
            for i2 = 1 : ngauss
                eta = nodes(i2);
                for i3 = 1 : ngauss
                    theta = nodes(i3);
                    [Floc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
                    Bloc = Bloc + Floc'*D*Floc*weights(i1)*weights(i2)*weights(i3)*detJ;
                    vol = vol + weights(i1)*weights(i2)*weights(i3)*detJ;
                end
            end
        end
        volumes(i) = vol;
        % Body force
        %for j = 1 : 8
        %    rhs0 = -vol/8*(divS(loc_coo(j,1),loc_coo(j,2),loc_coo(j,3)));
        %    ii = (3*(j-1)+1):(3*(j-1)+3);
        %    rhsloc(ii) = rhsloc(ii) + rhs0;
        %end

        % Global assembly
        loc_dof = 3*(loc_nod-1)+v3;
        loc_dof = loc_dof(:);
        [II,JJ] = meshgrid(loc_dof);
        Blist(k:k+575,:) = [JJ(:),II(:),Bloc(:)];
        k = k + 576;
    end

    B = sparse(Blist(:,1),Blist(:,2),Blist(:,3),3*nn,3*nn,size(Blist,1));
    B = 0.5*(B + B');

end
