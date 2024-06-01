function C = setStabilizationMat(C, nplas, tplas)

    ni = length(nplas);

    ID = 1:ni;

    % Open part does not require stabilization
    C(3*(ID(nplas)-1)+1,:) = 0;
    C(:,3*(ID(nplas)-1)+1) = 0;
    C(3*(ID(nplas)-1)+2,:) = 0;
    C(:,3*(ID(nplas)-1)+2) = 0;
    C(3*(ID(nplas)-1)+3,:) = 0;
    C(:,3*(ID(nplas)-1)+3) = 0;

    % Frictional components of the sliding part do not require stabilization
    C(3*(ID(tplas)-1)+2,:) = 0;
    C(:,3*(ID(tplas)-1)+2) = 0;
    C(3*(ID(tplas)-1)+3,:) = 0;
    C(:,3*(ID(tplas)-1)+3) = 0;

    % Save diagonal rows
    ID = 1:size(C,1);
    IDdiag = ID(sum(C~=0,1) == 1);
    Cdiag = C(IDdiag,IDdiag);

    % To force zero row sum
    C = C - diag(sum(C,1));

    % Restore diagonal rows
    C(IDdiag,IDdiag) = Cdiag;

end
