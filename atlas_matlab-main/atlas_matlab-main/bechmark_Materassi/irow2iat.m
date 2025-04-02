function iat = irow2iat(n,irow)

    nterm = length(irow);
    iat = zeros(n+1,1);
%     if (~issorted(irow))
%         irow = sort(irow);
%     end
%     if (irow(1) == 0)
%         irow = irow + 1;
%     end
    irow_old = 0;
    for k = 1 : nterm
        irow_new = irow(k);
        if (irow_new > irow_old)
            for j = irow_old+1 : irow_new
                iat(j) = k;
            end
            irow_old = irow_new;
        end
    end
    k = nterm + 1;
    for j = irow_old+1 : n+1
        iat(j) = k;
    end

end
