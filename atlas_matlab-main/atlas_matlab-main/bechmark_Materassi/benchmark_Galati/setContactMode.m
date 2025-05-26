function [F0,FRI] = setContactMode(FRIstick, FRIslide, F0stick, F0slide, masksP,nodePairsData)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Given the masks and the computed Jacobian and residual for all the
%     possible cases, sets the frictional part of the Jacobian to the
%     corresponding case.
%     In particular:
%         - if [Pt(u^{k+1})]_R+ < 0 --> Jacobian and residual = 0
%         - if the node is sliding --> take the sliding computation
%         - if the node is sticking --> take the sticking computation
%         - for the non-smooth cases:
%             - if []_R+ = 0 --> take 1/3 of slide+stick+0
%             - if []_(...) = 0 --> take 1/2 of slide+stick
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    v3 = [1;2;3];

    % Map interface nodes -> global dofs
    all_ntop = [nodePairsData.ntop]; % nodes at the interface

    nodesNneg = all_ntop(masksP.nneg);
    dofNneg = 3*(nodesNneg-1)+v3;
    dofNneg = dofNneg(:);
    ddofNneg = expand_dofs(dofNneg);

    nodesN0 = all_ntop(masksP.n0);
    dofN0 = 3*(nodesN0-1)+v3;
    dofN0 = dofN0(:);
    ddofN0 = expand_dofs(dofN0);

    nodesNpos = all_ntop(masksP.npos);
    dofNpos = 3*(nodesNpos-1)+v3;
    dofNpos = dofNpos(:);
    ddofNpos = expand_dofs(dofNpos);

    nodesTstick = all_ntop(masksP.tstick);
    dofTstick = 3*(nodesTstick-1)+v3;
    dofTstick = dofTstick(:);
    ddofTstick = expand_dofs(dofTstick);

    nodesTslide = all_ntop(masksP.tslide);
    dofTslide = 3*(nodesTslide-1)+v3;
    dofTslide = dofTslide(:);
    ddofTslide = expand_dofs(dofTslide);

    nodesT0 = all_ntop(masksP.t0);
    dofT0 = 3*(nodesT0-1)+v3;
    dofT0 = dofT0(:);
    ddofT0 = expand_dofs(dofT0);

    % Set mode for the Jacobian
    ndof = size(FRIstick,1);
    use_stick = true(ndof,1);
    use_slide = true(ndof,1);
    use_stick(ddofTslide) = false;
    use_slide(ddofTstick) = false;

    % Idea: sum the two cases (slide/stick), since:
    % - if [Pn]_+ < 0 --> it will be zeroed out
    % - if [Pn]_+ = 0 --> non-smooth case: sum 0 and divide by three
    % - if ||Pt|| = phi*Pn --> non-smooth case: just divide by 2
    FRI = sparse(ndof, ndof);
    FRI = FRI + FRIstick .* (use_stick * use_stick') + FRIslide .* (use_slide * use_slide');

    FRI(ddofNneg,ddofNneg) = 0.;
%     FRI(:,ddofNneg) = 0.;
%     FRI(ddofNneg,:) = 0.;
    
    % Handle the non-smooth cases
    FRI(ddofT0,ddofT0) = FRI(ddofT0,ddofT0)/2;
    FRI(ddofN0,ddofN0) = FRI(ddofN0,ddofN0)/3;
%     FRI(:,ddofT0) = FRI(:,ddofT0)/2;
%     FRI(:,ddofN0) = FRI(:,ddofN0)/3;
%     FRI(ddofT0,:) = FRI(ddofT0,:)/2;
%     FRI(ddofN0,:) = FRI(ddofN0,:)/3;

    % Set mode for the residual
    F0 = zeros(size(F0stick));
    F0(ddofN0) = 0.;
    F0(ddofNneg) = 0.;
    F0(ddofTstick) = F0stick(ddofTstick);
    F0(ddofTslide) = F0slide(ddofTslide);

end