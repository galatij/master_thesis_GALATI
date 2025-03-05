function [B,detJ] = cpt_shape(coord,csi,eta,theta,varargin)

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

    J = zeros(3,3);
    J(1,:) = dN*coord(:,1);
    J(2,:) = dN*coord(:,2);
    J(3,:) = dN*coord(:,3);
    detJ = abs(det(J));

    Btmp = J'\dN;

    B = zeros(6,24);
    for i = 1 : 8
        B(1,3*(i-1)+1) = Btmp(1,i);
        B(2,3*(i-1)+2) = Btmp(2,i);
        B(3,3*(i-1)+3) = Btmp(3,i);
        B(4,3*(i-1)+1) = Btmp(2,i);
        B(4,3*(i-1)+2) = Btmp(1,i);
        B(5,3*(i-1)+1) = Btmp(3,i);
        B(5,3*(i-1)+3) = Btmp(1,i);
        B(6,3*(i-1)+2) = Btmp(3,i);
        B(6,3*(i-1)+3) = Btmp(2,i);
    end
    
    if ~isempty(varargin)
        id = varargin{1};
        I = 1:3;
        dN = dN(I~=id,:);
        dN = dN*coord;
        detJ = norm(cross(J(1,:),J(2,:)));
    end
end
