function [M,ntot,rhs] = assemble_global(nn,nlam,K,Bt,rhs)

    ntot = 3*nn + 3*nlam;

    B = Bt';
    C = sparse(3*nlam,3*nlam);
    M = [K,B;Bt,C];
    rhs = [rhs;zeros(3*nlam,1)];

end
