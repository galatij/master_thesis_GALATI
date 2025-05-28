% function [F0,FRI] = setContactMode(FRIstick, FRIslide, F0stick, F0slide, masksP,nodePairsData)
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     Given the masks and the computed Jacobian and residual for all the
% %     possible cases, sets the frictional part of the Jacobian to the
% %     corresponding case.
% %     In particular:
% %         - if [Pt(u^{k+1})]_R+ < 0 --> Jacobian and residual = 0
% %         - if the node is sliding --> take the sliding computation
% %         - if the node is sticking --> take the sticking computation
% %         - for the non-smooth cases:
% %             - if []_R+ = 0 --> take 1/3 of slide+stick+0
% %             - if []_(...) = 0 --> take 1/2 of slide+stick
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Map interface nodes -> global dofs
%     all_ntop = [nodePairsData.ntop]; % nodes at the interface
% 
%     dofNneg = expand_dofs(all_ntop(masksP.nneg));
%     dofN0 = expand_dofs(all_ntop(masksP.n0));
%     dofNpos = expand_dofs(all_ntop(masksP.npos));
%     dofTstick = expand_dofs(all_ntop(masksP.tstick));
%     dofTslide = expand_dofs(all_ntop(masksP.tslide));
%     dofT0 = expand_dofs(all_ntop(masksP.t0));
% 
%     % Set mode for the Jacobian
%     ndof = size(FRIstick,1);
%     use_stick = true(ndof,1);
%     use_slide = true(ndof,1);
%     use_stick(dofTslide) = false;
%     use_slide(dofTstick) = false;
% 
%     % Idea: sum the two cases (slide/stick), since:
%     % - if [Pn]_+ < 0 --> it will be zeroed out
%     % - if [Pn]_+ = 0 --> non-smooth case: sum 0 and divide by three
%     % - if ||Pt|| = phi*Pn --> non-smooth case: just divide by 2
%     FRI = sparse(ndof, ndof);
%     FRI = FRI + FRIstick .* (use_stick * use_stick') + FRIslide .* (use_slide * use_slide');
% 
%     FRI(dofNneg,dofNneg) = 0.;
% %     FRI(:,dofNneg) = 0.;
% %     FRI(dofNneg,:) = 0.;
%     
%     % Handle the non-smooth cases
%     FRI(dofT0,dofT0) = FRI(dofT0,dofT0)/2;
%     FRI(dofN0,dofN0) = FRI(dofN0,dofN0)/3;
% %     FRI(:,dofT0) = FRI(:,dofT0)/2;
% %     FRI(:,dofN0) = FRI(:,dofN0)/3;
% %     FRI(dofT0,:) = FRI(dofT0,:)/2;
% %     FRI(dofN0,:) = FRI(dofN0,:)/3;
% 
%     % Set mode for the residual
%     F0 = zeros(size(F0stick));
%     F0(dofN0) = 0.;
%     F0(dofNneg) = 0.;
%     F0(dofTstick) = F0stick(dofTstick);
%     F0(dofTslide) = F0slide(dofTslide);
% 
% end


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

    % Map interface nodes -> global dofs
    all_ntop = [nodePairsData.ntop]; % nodes at the interface

    dofNneg = expand_dofs(all_ntop(masksP.nneg));
    dofN0 = expand_dofs(all_ntop(masksP.n0));
    dofNpos = expand_dofs(all_ntop(masksP.npos));
    dofTstick = expand_dofs(all_ntop(masksP.tstick));
    dofTslide = expand_dofs(all_ntop(masksP.tslide));
    dofT0 = expand_dofs(all_ntop(masksP.t0));

    % Set mode for the Jacobian
    ndof = size(FRIstick,1);
    use_stick = false(ndof,1);
    use_slide = false(ndof,1);
    non_smooth1 = false(ndof,1);
    non_smooth2 = false(ndof,1);
    use_stick(dofTstick) = true;
    use_slide(dofTslide) = true;
    non_smooth1(dofN0) = true;
    non_smooth2(dofT0) = true;

    FRI = sparse(ndof, ndof);
    FRI = FRI + FRIstick .* (use_stick * use_stick') ...
        + FRIslide .* (use_slide * use_slide') ...
        + (FRIstick+ FRIslide).*(non_smooth1 * non_smooth1')/3 ...
        + (FRIstick+ FRIslide).*(non_smooth2 * non_smooth2')/2;

    % Set mode for the residual
    F0 = zeros(size(F0stick));
    F0(dofN0) = 0.;
    F0(dofNneg) = 0.;
    F0(dofTstick) = F0stick(dofTstick);
    F0(dofTslide) = F0slide(dofTslide);

end