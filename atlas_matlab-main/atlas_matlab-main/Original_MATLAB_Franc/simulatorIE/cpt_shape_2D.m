function [N,detJ] = cpt_shape_2D(coord,csi,eta,theta,varargin)

    N = zeros(length(csi),8);
    N(:,1) = (1.0-csi).*(1.0-eta).*(1.0-theta)/8.0;
    N(:,2) = (1.0+csi).*(1.0-eta).*(1.0-theta)/8.0;
    N(:,3) = (1.0+csi).*(1.0+eta).*(1.0-theta)/8.0;
    N(:,4) = (1.0-csi).*(1.0+eta).*(1.0-theta)/8.0;
    N(:,5) = (1.0-csi).*(1.0-eta).*(1.0+theta)/8.0;
    N(:,6) = (1.0+csi).*(1.0-eta).*(1.0+theta)/8.0;
    N(:,7) = (1.0+csi).*(1.0+eta).*(1.0+theta)/8.0;
    N(:,8) = (1.0-csi).*(1.0+eta).*(1.0+theta)/8.0;

    if (nargout == 2)
        dN = zeros(3,8);
        dN(1,1) = -(1.0-eta)*(1.0-theta)/8.0;
        dN(1,2) = +(1.0-eta)*(1.0-theta)/8.0;
        dN(1,3) = +(1.0+eta)*(1.0-theta)/8.0;
        dN(1,4) = -(1.0+eta)*(1.0-theta)/8.0;
        dN(1,5) = -(1.0-eta)*(1.0+theta)/8.0;
        dN(1,6) = +(1.0-eta)*(1.0+theta)/8.0;
        dN(1,7) = +(1.0+eta)*(1.0+theta)/8.0;
        dN(1,8) = -(1.0+eta)*(1.0+theta)/8.0;
        dN(2,1) = -(1.0-csi)*(1.0-theta)/8.0;
        dN(2,2) = -(1.0+csi)*(1.0-theta)/8.0;
        dN(2,3) = +(1.0+csi)*(1.0-theta)/8.0;
        dN(2,4) = +(1.0-csi)*(1.0-theta)/8.0;
        dN(2,5) = -(1.0-csi)*(1.0+theta)/8.0;
        dN(2,6) = -(1.0+csi)*(1.0+theta)/8.0;
        dN(2,7) = +(1.0+csi)*(1.0+theta)/8.0;
        dN(2,8) = +(1.0-csi)*(1.0+theta)/8.0;
        dN(3,1) = -(1.0-csi)*(1.0-eta)/8.0;
        dN(3,2) = -(1.0+csi)*(1.0-eta)/8.0;
        dN(3,3) = -(1.0+csi)*(1.0+eta)/8.0;
        dN(3,4) = -(1.0-csi)*(1.0+eta)/8.0;
        dN(3,5) = +(1.0-csi)*(1.0-eta)/8.0;
        dN(3,6) = +(1.0+csi)*(1.0-eta)/8.0;
        dN(3,7) = +(1.0+csi)*(1.0+eta)/8.0;
        dN(3,8) = +(1.0-csi)*(1.0+eta)/8.0;

        id = varargin{1};
        I = 1:3;
        dN = dN(I~=id,:);
        J = dN*coord;
        detJ = norm(cross(J(1,:),J(2,:)));
    end

end
