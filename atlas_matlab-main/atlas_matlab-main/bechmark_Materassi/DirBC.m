function [A,rhs] = DirBC(ndir,inddir,presc,A,rhs)

    % Matlab sparse matrix are CSC

    [nrows,ncols] = size(A);
    [jcol,irow,coef] = find(A);                 % [irow,jcol,coef] = find(A);

    % Check if there are null diagonal terms and count them
    full_set = [1:nrows]';
    null_diag = setdiff(full_set,jcol(jcol==irow));
    nnull = length(null_diag);

    if (nnull > 0)
        % Update jcol, irow and coef
        i = [jcol;null_diag];               % [irow;null_diag];
        j = [irow;null_diag];               % [jcol;null_diag];
        aa = [coef;zeros(nnull,1)];

        % Save all arrays in a matrix
        P = [i,j,aa];

        % Sort new rows and columns
        P = sortrows(P,[2,1]);

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
        jj = inddir(i);     % index of the row where to apply b.c.
        WI(jj) = i;         % select the index of the b.c. from presc
    end

    for i = 1 : nrows
        if (WI(i) > 0)      % on row i I have to prescribe a Dirichlet b.c.
            k = WI(i);      % take the idx of the b.c. from presc
            j = iat(i);     % shift to CSR format
            while (jcol(j) < i)
                jcol_ind = jcol(j);
                rhs(jcol_ind) = rhs(jcol_ind)-coef(j)*presc(k);
                coef(j) = zero;
                j = j + 1;
            end
            % TODO: modify (not one)
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
