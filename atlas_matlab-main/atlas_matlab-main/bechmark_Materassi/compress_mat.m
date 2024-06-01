function A = compress_mat(A, blkSize)

    [iA, jA, cA] = find(A);
    iA = floor((iA-1)/blkSize) + 1;
    jA = floor((jA-1)/blkSize) + 1;
    oA = ones(nnz(A),1);
    A = sparse(iA, jA, cA);
    C = sparse(iA, jA, oA);
    [iA, jA, cA] = find(A);
    [~, ~, cC] = find(C);
    cA = cA./cC;
    A = sparse(iA, jA, cA);

end
