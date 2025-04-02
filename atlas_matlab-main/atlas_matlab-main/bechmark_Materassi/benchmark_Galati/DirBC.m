function [A,rhs] = DirBC(ndir,inddir,presc,A,rhs)

    % Matlab sparse matrix are CSC

    [nrows,ncols] = size(A);
    [jcol,irow,coef] = find(A);     %  why not [irow,jcol,coef]?

    % Check if there are null diagonal terms and count them
    full_set = [1:nrows]';
    null_diag = setdiff(full_set,jcol(jcol==irow));
    nnull = length(null_diag);

    if (nnull > 0)
        % Update jcol, irow and coef
        i = [jcol;null_diag];
        j = [irow;null_diag];
        aa = [coef;zeros(nnull,1)];

        % Save all arrays in a matrix
        P = [i,j,aa];

        % Sort new rows and columns
        P = sortrows(P,[2,1]);          % first sort rows by col (2), if equal sort by row (1)

        jcol = P(:,1);
        irow = P(:,2);
        coef = P(:,3);
    end

    % The CSR matrix is iat (size n+1), jcol (size nnz), coef (size nnz)
    iat = irow2iat(nrows,irow);

    WI = zeros(nrows,1);
    one = 1.0;
    zero = 0.0;

    for i = 1 : ndir
        jj = inddir(i);
        WI(jj) = i;
    end

    for i = 1 : nrows
        if (WI(i) > 0)
            k = WI(i);
            j = iat(i);
            while (jcol(j) < i)
                jcol_ind = jcol(j);
                rhs(jcol_ind) = rhs(jcol_ind)-coef(j)*presc(k);
                coef(j) = zero;
                j = j + 1;
            end
            coef(j) = one;
            for jj = j+1 : iat(i+1)-1
                jcol_ind = jcol(jj);
                rhs(jcol_ind) = rhs(jcol_ind)-coef(jj)*presc(k);
                coef(jj) = zero;
            end
        else
            for j = iat(i) : iat(i+1)-1
                jcol_ind = jcol(j);
                if (WI(jcol_ind) > 0)
                    coef(j) = zero;
                end
            end
        end
    end

    for i = 1 : ndir
        jj = inddir(i);
        rhs(jj) = presc(i);
    end

    % Assemble the final matrix
    A = sparse(jcol,irow,coef,nrows,ncols);

end
