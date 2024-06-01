function [int_area] = cpt_area_int(ngauss, coord, topol, elem, list)

    loc_coo = coord(topol(elem,:),:);

    csi   = [-1;+1;+1;-1;-1;+1;+1;-1];
    eta   = [-1;-1;+1;+1;-1;-1;+1;+1];
    theta = [-1;-1;-1;-1;+1;+1;+1;+1];

    xi = [csi,eta,theta];
    Nloc = cpt_shape_2D(loc_coo,csi,eta,theta);
    X = ismember(Nloc*loc_coo,coord(list,:),'row');
    id = 0;
    for i = 1 : 3
        if (std(xi(X,i)) == 0)
            xi_id = i;
            xi_val = mean(xi(X,i));
        end
    end

    [nodes,weights] = gausspoints(ngauss);

    ID0 = [1:3]';
    ID = zeros(3,1);
    ID(ID0~=xi_id) = [1:2];

    int_area = zeros(1,4);
    tmp = zeros(3,1);
    tmp(xi_id) = xi_val;
    for i1 = 1 : ngauss
        csi = nodes(i1);
        tmp(find(ID==1)) = csi;
        for i2 = 1 : ngauss
            eta = nodes(i2);
            tmp(find(ID==2)) = eta;
            [Nloc,detJ] = cpt_shape_2D(loc_coo,tmp(1),tmp(2),tmp(3),xi_id);
            int_area = int_area + weights(i1)*weights(i2)*detJ*Nloc(X);
        end
    end

end
