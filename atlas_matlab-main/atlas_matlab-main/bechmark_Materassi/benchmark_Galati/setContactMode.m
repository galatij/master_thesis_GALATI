function [F0,FRI] = setContactMode(FRIstick, FRIslide, F0stick, F0slide, masksP,nodePairsData)

    v3 = [1;2;3];

    % map nni -> global dofs
    all_ntop = [nodePairsData.ntop];

    nodesNneg = all_ntop(masksP.nneg);
    dofNneg = 3*(nodesNneg-1)+v3;
    dofNneg = dofNneg(:);

    nodesN0 = all_ntop(masksP.n0);
    dofN0 = 3*(nodesN0-1)+v3;
    dofN0 = dofN0(:);

    nodesNpos = all_ntop(masksP.npos);
    dofNpos = 3*(nodesNpos-1)+v3;
    dofNpos = dofNpos(:);

    nodesTstick = all_ntop(masksP.tstick);
    dofTstick = 3*(nodesTstick-1)+v3;
    dofTstick = dofTstick(:);

    nodesTslide = all_ntop(masksP.tslide);
    dofTslide = 3*(nodesTslide-1)+v3;
    dofTslide = dofTslide(:);

    nodesT0 = all_ntop(masksP.t0);
    dofT0 = 3*(nodesT0-1)+v3;
    dofT0 = dofT0(:);

    % Set mode for the Jacobian
    ndof = size(FRIstick,1);
    use_stick = true(ndof,1);
    use_slide = true(ndof,1);
    use_stick(dofTslide) = false;
    use_slide(dofTstick) = false;

    FRI = sparse(ndof, ndof);
    FRI = FRI + FRIstick .* (use_stick * use_stick') + FRIslide .* (use_slide * use_slide');

    FRI(dofNneg,dofNneg) = 0.;

    % Handle the non-smooth cases
    FRI(dofT0,dofT0) = FRI(dofT0,dofT0)/2;
    FRI(dofN0,dofN0) = FRI(dofN0,dofN0)/3;

    % Set mode for the residual
    F0 = zeros(size(F0stick));
    F0(dofN0) = 0.;
    F0(dofNneg) = 0.;
    F0(dofTstick) = F0stick(dofTstick);
    F0(dofTslide) = F0slide(dofTslide);

end