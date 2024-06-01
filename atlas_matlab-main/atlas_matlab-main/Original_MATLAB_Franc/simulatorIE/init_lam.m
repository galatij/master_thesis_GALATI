function [lam_ini] = init_lam(nlam,lam0)

    lam_ini = zeros(3*nlam,1);
    lam_ini(1:3:3*nlam) = lam0;

end
