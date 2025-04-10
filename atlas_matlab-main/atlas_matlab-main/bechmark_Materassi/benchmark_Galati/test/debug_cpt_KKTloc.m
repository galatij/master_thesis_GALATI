clear; clc;
addpath("../")

E = 1;
nu = 0.3;
D = cpt_elas_mat(E, nu);
ngauss = 2;
gamma = 10;
elem = 1;
alpha = 1;
% Simple hexahedral element
coord = [
    0, 0, 0; 
    1, 0, 0;
    1, 1, 0;
    0, 1, 0;
    0, 0, 1;
    1, 0, 1;
    1, 1, 1;
    0, 1, 1
];

topol = [1, 2, 3, 4, 5, 6, 7, 8];

% Ass: contact at x = 0
face_nodes = [1, 4, 8, 5];
normal = [-1; 0; 0];
S_n = normal'*cpt_normal(normal);
nN = repmat(normal',1,8);

%% From cpt_KKTloc.m
csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
theta = [-1;-1;-1;-1;+1;+1;+1;+1];
    
% Identify the fixed direction
xi = [csi,eta,theta];
Nloc_top = cpt_shape_2D(coord,csi,eta,theta);
X_top = ismember(Nloc_top*coord,coord(face_nodes,:),'row');
id = 0;

for i = 1 : 3
    if (std(xi(X_top,i)) == 0)
        xi_id_top = i;
        xi_val_top = mean(xi(X_top,i));
    end
end

[nodes,weights] = gausspoints(ngauss);


ID0 = (1:3)';
ID_top = zeros(3,1);
ID_top(ID0~=xi_id_top) = 1:2;
    
    
% Compute local contribuition on the top face (biased formulation)
KKTloc = zeros(24,24);
tmp_top = zeros(3,1);
tmp_top(xi_id_top) = xi_val_top;
for i1 = 1 : ngauss
    csi = nodes(i1);
    tmp_top(ID_top==1) = csi;
    for i2 = 1 : ngauss
        eta = nodes(i2);
        tmp_top(ID_top==2) = eta;
        [Bloc_top,detJ] = cpt_shape(coord,tmp_top(1),tmp_top(2),tmp_top(3),xi_id_top);        % shape derivatives 6x24
        [Nloc_top] = cpt_shape_2D(coord,tmp_top(1),tmp_top(2),tmp_top(3));                % shape functions 1x8
        Nloc_top = repelem(Nloc_top,1,3);                                      % shape functions 1x24
        disp('Shape functions (top face):');
        disp(Nloc_top);
                    
        P = gamma*Nloc_top.*nN - (S_n*D*Bloc_top);                 % P 1x24
        Palpha = gamma*Nloc_top.*nN - alpha*(S_n*D*Bloc_top);
        
        sol_test = coord';
        sol_test = sol_test(:);
        disp(Nloc_top(:).*sol_test.*nN');
        disp('P*coord:')
        disp(P*sol_test);      % works correctly
    
        % KKTloc multiplies the top displacement
        KKTloc = KKTloc + 1/gamma*P'*Palpha*weights(i1)*weights(i2)*detJ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computed local KKTloc:')
disp(KKTloc)
spy(KKTloc)

% Symmetry check
if norm(KKTloc - KKTloc', 'fro') < 1e-12
    disp('KKTloc is symmetric');
else
    warning('KKTloc is NOT symmetric');
end
