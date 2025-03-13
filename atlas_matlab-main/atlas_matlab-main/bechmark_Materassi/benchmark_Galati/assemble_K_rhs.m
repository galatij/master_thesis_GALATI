function [K,rhs,E,volumes] = assemble_K_rhs(ngauss,nn,coord,ne,topol,E0,nu,e2f,faces,faceData,bound)
                                  
    [nodes,weights] = gausspoints(ngauss);

    if (length(E0) == 1)
        E = E0*ones(ne,1);
    else
        E = E0;
    end

    f2e = e2f';

    k = 1;
    Klist = zeros(ne*576,3);
    rhs = zeros(3*nn,1);
    v3 = [1;2;3];
    face = [1,2,3,4;1,2,6,5;2,3,7,6;4,3,7,8;1,4,8,5;5,6,7,8];
    volumes = zeros(ne,1);

    for i = 1 : ne
        % Local stiffness matrix
        D = cpt_elas_mat(E(i),nu);

        loc_nod = topol(i,:);
        loc_coo = coord(loc_nod,:);
        Kloc = zeros(24,24);
        rhsloc = zeros(24,1);
        vol = 0.0;
        for i1 = 1 : ngauss
            csi = nodes(i1);
            for i2 = 1 : ngauss
                eta = nodes(i2);
                for i3 = 1 : ngauss
                    theta = nodes(i3);
                    [Bloc,detJ] = cpt_shape(loc_coo,csi,eta,theta);
                    Kloc = Kloc + Bloc'*D*Bloc*weights(i1)*weights(i2)*weights(i3)*detJ;
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


        % Check if the element is at the interface
%         loc_interf = find(interf2e(:,i));
%         if ~isempty(loc_interf)
%             %% Check if the element is a top or bottom element
%             if (interfData(loc_interf).etop == i)
%                 list = interfData(loc_interf).top;
%                 normal = interfData(loc_interf).normal;
%             else
%                 list = interfData(loc_interf).bottom;
%                 normal = interfData(loc_interf).normal;
%             end
%             N = normal*normal';
%             Gloc = cpt_shape_interf(ngauss, coord, topol, i, list, ...
%                                    N, gamma, D);
%             Kloc = Kloc + Gloc;
%         end
        

        % Neumann boundary condition
        loc_faces = find(f2e(:,i));
        for jf = 1 : length(bound)
            for kf = 1 : 6
                if (bound(jf).isbound(loc_faces(kf)) == 1)
                    loc_nodes = faces(loc_faces(kf),:);
                    iif = find(loc_faces(kf)==bound(jf).ID);
                    val = bound(jf).values(iif,:);
                    val = val(:);
                    area = faceData(i).area(kf);
                    int_area = cpt_area_int(ngauss, coord, topol, i, loc_nodes);
                    jjf = 3*(loc_nodes-1)+v3;            % matrix
                    jjf = jjf(:);                        % hstack
                    rhs0 = val*int_area;
                    rhs(jjf) = rhs(jjf) + rhs0(:);
                end
            end
        end

        % Global assembly
        loc_dof = 3*(loc_nod-1)+v3;
        loc_dof = loc_dof(:);
        [II,JJ] = meshgrid(loc_dof);
        Klist(k:k+575,:) = [JJ(:),II(:),Kloc(:)];
        k = k + 576;
        rhs(loc_dof) = rhs(loc_dof) + rhsloc;
    end

    K = sparse(Klist(:,1),Klist(:,2),Klist(:,3),3*nn,3*nn,size(Klist,1));
    K = 0.5*(K + K');

end
