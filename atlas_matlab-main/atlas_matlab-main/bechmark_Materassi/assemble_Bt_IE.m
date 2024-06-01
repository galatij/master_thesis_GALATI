function [Bt] = assemble_Bt_IE(ngauss,coord,topol,interf,interfData)

    % With numerical integration of the basis functions restricted on the face

    ni = size(interf,1);
    nn = size(coord,1);
    irow = zeros(72*ni,1);
    jcol = zeros(72*ni,1);
    coef = zeros(72*ni,1);
    v3 = [1;2;3];
    k = 1;
    for i = 1 : ni

        R = interfData(i).R;
        int_area = cpt_area_int(ngauss, coord, topol, interfData(i).etop, interfData(i).top);

        list = interfData(i).top;
        for j = 1 : 4
            irow(k:k+2) = 3*(i-1)+1;
            jcol(k:k+2) = 3*(list(j)-1)+v3;
            irow(k+3:k+5) = 3*(i-1)+2;
            jcol(k+3:k+5) = 3*(list(j)-1)+v3;
            irow(k+6:k+8) = 3*(i-1)+3;
            jcol(k+6:k+8) = 3*(list(j)-1)+v3;
            coef(k:k+8) = -int_area(j)*R(:);
            k = k + 9;
        end

        list = interfData(i).bottom;
        for j = 1 : 4
            irow(k:k+2) = 3*(i-1)+1;
            jcol(k:k+2) = 3*(list(j)-1)+v3;
            irow(k+3:k+5) = 3*(i-1)+2;
            jcol(k+3:k+5) = 3*(list(j)-1)+v3;
            irow(k+6:k+8) = 3*(i-1)+3;
            jcol(k+6:k+8) = 3*(list(j)-1)+v3;
            coef(k:k+8) = int_area(j)*R(:);
            k = k + 9;
        end
    end
    Bt = sparse(irow,jcol,coef,3*ni,3*nn);

end
