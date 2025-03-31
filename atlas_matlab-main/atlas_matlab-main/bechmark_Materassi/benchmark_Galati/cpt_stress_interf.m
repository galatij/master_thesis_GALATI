function [sigma_n, sigma_t] = cpt_stress_interf(stress, nodePairsData)
    
    nni = numel(nodePairsData);
    sigma_n = zeros(nni,1);
    sigma_t = zeros(nni,2);
    for i = 1:nni
        ntop = nodePairsData(i).ntop;
        loc_stress = stress(ntop,:)';
        n = nodePairsData(i).normal;
        N = cpt_normal(n);

        % Compute two orthonormal tangential basis vectors
        if abs(n(2)) < 1e-6 && abs(n(3)) < 1e-6  % If n is aligned with z-axis
            e1 = [0;1;0];  % Pick a different arbitrary vector
        else
            e1 = [1;0;0];  
        end
        t1 = cross(e1, n);
        t1 = t1 / norm(t1);
        t2 = cross(n, t1);
        t2 = t2 / norm(t2);

        sigma_n(i) = n'*N*loc_stress;
        sigma_t(i,1) = t1'*N*loc_stress;
        sigma_t(i,2) = t2'*N*loc_stress;
    end

end


% function [sigma_n, sigma_t] = cpt_stress_interf(nodePairsData,sol)
%     
%     nni = numel(nodePairsData);
%     v3 = [1;2;3];
%     stress = zeros(nni,6);
%     sigma_n = zeros(nni,1);
%     sigma_t = zeros(nni,2);
%     
%     for i = 1 : nni
%         % Just compute the stress on the top face
%         D = nodePairsData(i).Dtop;
%         loc_coo = nodePairsData(i).coord;
%         loc_nod = nodePairsData(i).ntop;
%         loc_dof = 3*(loc_nod-1)+v3;
%         loc_dof = loc_dof(:);
%         loc_sol = sol(loc_dof);
% 
%         % Compute geometric vectors for normal and tangential
%         % components
%         n = nodePairsData(i).normal;
%         N = cpt_normal(n);
% 
%         % TODO: generalize for general n
%         % Compute two orthonormal tangential basis vectors
%         if abs(n(1)) < 1e-6 && abs(n(2)) < 1e-6  % If n is aligned with z-axis
%             e1 = [0;1;0];  % Pick a different arbitrary vector
%         else
%             e1 = [1;0;0];  
%         end
%         t1 = cross(e1, n);
%         t1 = t1 / norm(t1);
%         t2 = cross(n, t1);
%         t2 = t2 / norm(t2);
%         
%         % Compute the stress
%         Bloc = cpt_shape(loc_coo,0,0,0);
%         loc_eps = Bloc*loc_sol;             % [exx; eyy; ezz; gammaxz; gammaxz; gammayz]
%         loc_stress = D*loc_eps;             % [sxx; syy; szz; sxz; sxz; syz]
%         stress(i,:) = loc_stress;
%         sigma_n(i) = n'*N*loc_stress;
%         sigma_t(i,1) = t1'*N*loc_stress;
%         sigma_t(i,2) = t2'*N*loc_stress;
%         % TODO: check if to do the same for the bottom face (unbiased
%         % formulation)
%     end
% end