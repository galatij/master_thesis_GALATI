function [D] = cpt_elas_mat(E,nu)

    % Constituent matrix
    D = zeros(6,6);
    D(1,1) = 1-nu;
    D(1,2) = nu;
    D(1,3) = nu;
    D(2,1) = nu;
    D(2,2) = 1-nu;
    D(2,3) = nu;
    D(3,1) = nu;
    D(3,2) = nu;
    D(3,3) = 1-nu;
    D(4,4) = (1-2*nu)/2;
    D(5,5) = (1-2*nu)/2;
    D(6,6) = (1-2*nu)/2;
    D = E/((1+nu)*(1-2*nu))*D;

end
