function C = setStabilizationMat(C, nplas, tplas)

    ni = length(nplas);

    ID = 1:ni;
    ID_nplas = ID(nplas);
    ID_tplas = ID(tplas);

    % Open part and sliding parts do not require stabilization
    C(:,3*(ID_nplas-1)+1) = 0;
    C(:,3*(ID_tplas-1)+2) = 0;
    C(:,3*(ID_tplas-1)+3) = 0;
    C = C';
    C(:,3*(ID_nplas-1)+1) = 0;
    C(:,3*(ID_tplas-1)+2) = 0;
    C(:,3*(ID_tplas-1)+3) = 0;

    % Get the IDs of the rows having only the diagonal entry
    nnzr = sum(spones(C));
    ind_diag = find(nnzr==1);
    % Save those diagonal entries
    diagC = full(diag(C));
    D = zeros(size(C,1),1);
    D(ind_diag) = diagC(ind_diag);

    % Force zero row sum
    C = C - diag( sum(C)' - sparse(D) );

end
