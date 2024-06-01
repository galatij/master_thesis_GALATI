function Bt = setCouplingMat(Bt, nplas, tplas)

    ni = length(nplas);

    ID = 1:ni;

    Bt(3*(ID(nplas)-1)+1,:) = 0;
    Bt(3*(ID(tplas)-1)+2,:) = 0;
    Bt(3*(ID(tplas)-1)+3,:) = 0;

end
