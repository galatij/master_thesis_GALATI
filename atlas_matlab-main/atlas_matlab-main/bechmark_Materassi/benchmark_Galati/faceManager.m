function [nf, faces, faceData, e2f, areaf, interfData] = ...
    faceManager(ne, topol, ni, interf, coord, interf2e, bar_ele)

    e2i = interf2e';
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
            faces(6*(i-1)+j,:) = topol(i,face(j,:));
            facesS(6*(i-1)+j,:) = sort(topol(i,face(j,:)));
            k = k + 1;
            irow(k) = i;
            jcol(k) = 6*(i-1)+j;
        end
    end
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
        elem = find(e2f(:,i));
        % Loop over elements sharing face i
        for j = 1 : length(elem)
            % Get list of faces for each element
            jj = find(f2e(:,elem(j)));
            % Find index of face i in the face list of each element
            kk = find(jj==i);
            % Get element barycentre
            x0 = bar_ele(elem(j),:)';
            % Check normal orientation
            if ((loc_bar-x0)'*normf0 > 0)
                faceData(elem(j)).normal(kk,:) = normf0(:);
            else
                faceData(elem(j)).normal(kk,:) = -normf0(:);
            end
            faceData(elem(j)).area(kk) = areaf(i);
            faceData(elem(j)).list = jj(:);
        end
    end
    
    interfData = repmat(struct('top', zeros(4,1), 'bottom', zeros(4,1), 'etop', 1, 'ebottom', 1, ...
        'area', 1, 'normal', zeros(3,1), 'bar', zeros(3,1), 'normEdges', zeros(4,3)), ni, 1);
    for i = 1 : ni
        iloc = interf(i,1:4);
        v1 = coordT(:,iloc(2)) - coordT(:,iloc(1));
        v2 = coordT(:,iloc(3)) - coordT(:,iloc(1));
        nrm1 = cross(v1, v2);
        area1 = norm(nrm1);
        v1 = coordT(:,iloc(4)) - coordT(:,iloc(3));
        v2 = coordT(:,iloc(1)) - coordT(:,iloc(3));
        nrm2 = cross(v1, v2);
        area2 = norm(nrm2);
        interfData(i).top = interf(i,1:4);
        interfData(i).bottom = interf(i,5:8);
        interfData(i).area = 0.5*(area1 + area2);
        % Compute the normal
        normf0 = (nrm1 + nrm2) / (area1 + area2);
        % Get elements sharing a face with interf element i
        [ii, ~, dd] = find(e2i(:,i));
        itop = find(dd==1);
        ibot = find(dd==2);
        % Select the top element
        ielem = ii(itop);
        interfData(i).etop = ielem;
        interfData(i).ebottom = ii(ibot);
        loc_bar = mean(coordT(:,interf(i,1:4)),2);
        x0 = bar_ele(ielem,:)';
        % To ensure TOP -> BOTTOM normal orientation
        if ((loc_bar-x0)'*normf0 > 0)
            interfData(i).normal = normf0(:);
        else
            interfData(i).normal = -normf0(:);
        end
        % Compute local rotation matrix
        nrm = interfData(i).normal;
        [Q, ~] = qr(nrm);
        D = zeros(3,1);
        D(1) = nrm'*Q(:,1);
        for j = 2 : 3
            [~, pos] = max(abs(Q(:,j)));
            D(j) = sign(Q(pos,j));
        end
        interfData(i).R = Q * diag(D);
        % Save the barycenter
        interfData(i).bar = loc_bar;
      
    end

end
