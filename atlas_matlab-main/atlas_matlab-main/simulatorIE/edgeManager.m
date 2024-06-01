function [edges,edgeData,f2e,interfData] = edgeManager(interf,coord,interfData)
%-----------------------------------------------------------------------------------------
%
% Function to compute interface edge data
% edges(n_edge):        list of IE edges
% edgeData(:,n_edge):   struct for the edge data
%                           n1: node 1
%                           n2: node 2
%                           len: length
%                           normal: normal vector to the edge
%                           bar: baricenter of the edge
% f2e:                  mapping from faces to edge
% interfData:           this struct is only updated in the normEdges part. nromEdges are
%                       normals to the face edges pointing outward
%
%-----------------------------------------------------------------------------------------

   ni = size(interf,1);

   % face-edge topology
   edges = zeros(4*ni,2);
   irow = zeros(4*ni,1);
   jcol = zeros(4*ni,1);
   coef = ones(4*ni,1);
   k = 0;
   for i = 1 : ni
       for j = 1 : 4
           i1 = mod(j-1,4)+1;
           i2 = mod(j,4)+1;
           edges(4*(i-1)+j,:) = sort([interf(i,i1),interf(i,i2)]);
           k = k + 1;
           irow(k) = i;
           jcol(k) = 4*(i-1)+j;
       end
   end
   % Eliminate redundant edges
   [edges,~,idedges] = unique(edges,'rows');
   nedges = size(edges,1);
   jcol = idedges(jcol);
   % Create mapping
   f2e = sparse(irow,jcol,coef);
   e2f = f2e';

   % Create database of edge properties
   edgeData = repmat(struct( 'n1', 1, 'n2', 1, 'len', 1, 'normal', zeros(3,1),...
                             'bar', zeros(3,1) ),...
                             nedges, 1);

   coordEdge = zeros(nedges,3);
   for i = 1 : nedges
       coord1 = coord(edges(i,1),:);
       coord2 = coord(edges(i,2),:);
       edgeData(i).n1 = edges(i,1);
       edgeData(i).n2 = edges(i,2);
       edgeData(i).bar = (coord1 + coord2)/2.0;
       edgeData(i).len = norm(coord1 - coord2);
       % List of faces sharing edge i
       faces = find(f2e(:,i));
       if (length(faces) == 2)
           f1 = faces(1);
           f2 = faces(2);
           nrm1 = interfData(f1).normal;
           nrm2 = interfData(f2).normal;
           A1 = interfData(f1).area;
           A2 = interfData(f2).area;
           nrmf = (nrm1*A1 + nrm2*A2)/(A1+A2);
           nrme = cross(nrmf, coord1-coord2);
           nrme = nrme/norm(nrme);
       else
           f1 = faces(1);
           nrmf = interfData(f1).normal;
           nrme = cross(nrmf, coord1-coord2);
           nrme = nrme/norm(nrme);
       end
       edgeData(i).normal = nrme;

       loc_bar = edgeData(i).bar;
       loc_bar = loc_bar(:);
       % Loop over elements sharing face i
       for j = 1 : length(faces)
           % Get list of edges for each face
           jj = find(e2f(:,faces(j)));
           % Find index of edge i in the edge list of each face
           kk = find(jj==i);
           % Get face barycentre
           x0 = interfData(faces(j)).bar;
           % Check normal orientation (the edge normal points outward)
           if ((loc_bar-x0)'*nrme' > 0)
               interfData(faces(j)).normEdges(kk,:) = nrme(:);
           else
               interfData(faces(j)).normEdges(kk,:) = -nrme(:);
           end
       end
   end

end
