function [faces,faceData,e2f,areaf,interfData] = ...
                faceManager(elem,interf,coord,i2e,bar_ele);
%-----------------------------------------------------------------------------------------
%
% Creates face information
% faces(nf):          list of all the face in the mesh
% faceData(:,ne):     struct for the elements in the mesh
%                     list: list of faces for this element
%                     area: area of the faces for this element
%                     normal: normal vector of the faces for this element (pointing outside)
% e2f:                sparse matrix mapping elements to its faces
% areaf(nf):          area of the faces (redundant)
% interfData(:,ni)    struct for the IE in the mesh
%                     top: list of nodes of top
%                     bottom: list of nodes of bottom
%                     etop: element sharing top face
%                     ebotto: element sharing bottom face
%                     area: area of the IE
%                     normal: normal versor from BOTTOM to TOP
%                     bar: baricenter of the IE
%                     normEdges: normal to the edges pointing outside (loaded later)
%                     R_l2g: rotation matrix from local to global reference
%
%-----------------------------------------------------------------------------------------

    % Tolerance to check X-coaxiality
    tol_Xax = 1.e-9;

    ne = size(elem,1);
    ni = size(interf,1);
    e2i = i2e';
    coordT = coord';
    faces = zeros(6*ne,4);
    facesS = zeros(6*ne,4);
    irow = zeros(6*ne,1);
    jcol = zeros(6*ne,1);
    coef = ones(6*ne,1);
    face = [1,2,3,4;1,2,6,5;2,3,7,6;4,3,7,8;1,4,8,5;5,6,7,8];
    k = 0;
    for i = 1 : ne
        for j = 1 : 6
            faces(6*(i-1)+j,:) = elem(i,face(j,:));
            facesS(6*(i-1)+j,:) = sort(elem(i,face(j,:)));
            k = k + 1;
            irow(k) = i;
            jcol(k) = 6*(i-1)+j;
        end
    end
    % Eliminate redundant faces
    [~, idunique, idfaces] = unique(facesS, 'rows');
    faces = faces(idunique,:);
    jcol = idfaces(jcol);
    e2f = sparse(irow, jcol, coef);
    f2e = e2f';
    nf = size(faces,1);
    areaf = zeros(nf,1);
    faceData = repmat(struct('list', 1, 'area', zeros(6,1), 'normal', zeros(6,3)), ne, 1);
    for i = 1 : nf
        iloc = faces(i,:);
        v1 = coordT(:,iloc(2)) - coordT(:,iloc(1));
        v2 = coordT(:,iloc(3)) - coordT(:,iloc(1));
        nrm1 = cross(v1, v2);
        area1 = norm(nrm1);
        v1 = coordT(:,iloc(4)) - coordT(:,iloc(3));
        v2 = coordT(:,iloc(1)) - coordT(:,iloc(3));
        nrm2 = cross(v1, v2);
        area2 = norm(nrm2);
        areaf(i) = 0.5*(area1 + area2);
        normf0 = (nrm1 + nrm2) / (area1 + area2);
        loc_bar = mean(coordT(:,faces(i,:)),2);
        % List of elements sharing face i
        el_list = find(e2f(:,i));
        % Loop over elements sharing face i
        for j = 1 : length(el_list)
            % Get list of faces for each element
            jj = find(f2e(:,el_list(j)));
            % Find index of face i in the face list of each element
            kk = find(jj==i);
            % Get element barycentre
            x0 = bar_ele(el_list(j),:)';
            % Check normal orientation
            if ((loc_bar-x0)'*normf0 > 0)
                faceData(el_list(j)).normal(kk,:) = normf0(:);
            else
                faceData(el_list(j)).normal(kk,:) = -normf0(:);
            end
            faceData(el_list(j)).area(kk) = areaf(i);
            faceData(el_list(j)).list = jj(:);
        end
    end

    interfData = repmat( struct('top', zeros(4,1), 'bottom', zeros(4,1), 'etop', 1,...
                                'ebottom', 1, 'area', 1, 'normal', zeros(3,1),...
                                'bar', zeros(3,1), 'normEdges', zeros(4,3),...
                                'R_l2g' , zeros(3)),...
                                 ni , 1);
    for i = 1 : ni
        % BOTTOM
        iloc = interf(i,1:4);
        v1 = coordT(:,iloc(2)) - coordT(:,iloc(1));
        v2 = coordT(:,iloc(3)) - coordT(:,iloc(1));
        nrm1 = cross(v1, v2);
        area1 = norm(nrm1);
        v1 = coordT(:,iloc(4)) - coordT(:,iloc(3));
        v2 = coordT(:,iloc(1)) - coordT(:,iloc(3));
        nrm2 = cross(v1, v2);
        area2 = norm(nrm2);
        interfData(i).bottom = interf(i,1:4);
        interfData(i).top = interf(i,5:8);
        interfData(i).area = 0.5*(area1 + area2);
        % Compute the normal
        normf0 = (nrm1 + nrm2) / (area1 + area2);
        % Get elements sharing a face with interf element i
        [ii, ~, dd] = find(e2i(:,i));
        ibot = find(dd==1);
        itop = find(dd==2);
        % Select the top element
        ielem = ii(itop);
        interfData(i).etop = ielem;
        interfData(i).ebottom = ii(ibot);
        loc_bar = mean(coordT(:,interf(i,1:4)),2);
        x0 = bar_ele(ielem,:)'; % Baricenter of top element
        % To ensure BOTTOM -> TOP normal orientation
        if ((x0-loc_bar)'*normf0 > 0)
            interfData(i).normal = normf0(:);
        else
            interfData(i).normal = -normf0(:);
        end
        % Compute local rotation matrix
        N_v = interfData(i).normal;
        if N_v(2) < tol_Xax && N_v(3) < tol_Xax
           % The normal is aligned with X axis
           if N_v(1) > 0
              % Same direction
              interfData(i).R_l2g = eye(3);
           else
              % Opposite direction
              interfData(i).R_l2g = -eye(3);
              interfData(i).R_l2g(3,3) = 1.0;
           end
        else
           % The normal forms a non trivial angol with X axis
           T1 = cross(N_v,[1;0;0]);
           T2 = cross(N_v,T1);
           interfData(i).R_l2g = [N_v T1 T2];
        end
        % Save the barycenter
        interfData(i).bar = loc_bar;
    end

end
