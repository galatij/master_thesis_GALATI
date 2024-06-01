function x = direct_solver(A, b)

    x = A\b;
    return

    nA = size(A,1);
    D = full(diag(A));
    D(D==0) = 1;
    D = 1./sqrt(abs(D));
    D = spdiags(D,0,nA,nA);
    x = D*((D*A*D)\(D*b));
    norm(b-A*x)/norm(b)

end
