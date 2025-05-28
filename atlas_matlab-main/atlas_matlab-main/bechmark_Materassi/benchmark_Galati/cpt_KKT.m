function [C_k, KKT] = cpt_KKT(ngauss,coord,topol,E, nu, ...
                          interfData, nodePairsData, gamma, alpha, ...
                          dsol, masksP)
    
    TEST = true;
    ni = size(interfData,1);
    nn = size(coord,1);
    
    %% Evaluate the modulus
    v3 = [1;2;3];

    KKTlist11 = zeros(ni*576, 3);
    KKTlist12 = zeros(ni*576, 3);
    KKTlist21 = zeros(ni*576, 3);
    KKTlist22 = zeros(ni*576, 3);
    k = 1;
    for i = 1 : ni 

        % Compute contribution to the top face only (biased formulation)
        [KKT11,KKT12,KKT21,KKT22] = cpt_KKTloc(ngauss, coord, topol, interfData, ...
            i, E, nu, gamma, alpha, masksP);
        
        top_nod = topol(interfData(i).etop,:);
        bot_nod = topol(interfData(i).ebottom,:);
        top_dof = 3*(top_nod-1)+v3;
        top_dof = top_dof(:);
        bot_dof = 3*(bot_nod-1)+v3;
        bot_dof = bot_dof(:);

        [II_top,JJ_top] = meshgrid(top_dof);
        KKTlist11(k:k+575,:) = [JJ_top(:),II_top(:),KKT11(:)];
%         KKTlist11(k:k+575,:) = [II_top(:),JJ_top(:),KKT11(:)];

        [II_bot,JJ_top] = meshgrid(bot_dof, top_dof);
        KKTlist12(k:k+575,:) = [JJ_top(:),II_bot(:),KKT12(:)];
%         KKTlist12(k:k+575,:) = [II_bot(:),JJ_top(:),KKT12(:)];

        [II_top,JJ_bot] = meshgrid(top_dof, bot_dof);
        KKTlist21(k:k+575,:) = [JJ_bot(:),II_top(:),KKT21(:)];
%         KKTlist21(k:k+575,:) = [II_top(:),JJ_bot(:),KKT21(:)];

        [II_bot,JJ_bot] = meshgrid(bot_dof);
%         KKTlist22(k:k+575,:) = [JJ_bot(:),II_bot(:),KKT22(:)];
        KKTlist22(k:k+575,:) = [II_bot(:),JJ_bot(:),KKT22(:)];


        k = k + 576;
                
    end
    % CHECK: ...(:,1), ...(:,2) or the reverse?
    KKT = sparse(KKTlist11(:,2),KKTlist11(:,1),KKTlist11(:,3),3*nn,3*nn,size(KKTlist11,1)) ...
        + sparse(KKTlist12(:,2),KKTlist12(:,1),KKTlist12(:,3),3*nn,3*nn,size(KKTlist12,1)) ...
        + sparse(KKTlist21(:,2),KKTlist21(:,1),KKTlist21(:,3),3*nn,3*nn,size(KKTlist21,1)) ...
        + sparse(KKTlist22(:,2),KKTlist22(:,1),KKTlist22(:,3),3*nn,3*nn,size(KKTlist22,1));

    C_k = KKT * dsol;
    

end
