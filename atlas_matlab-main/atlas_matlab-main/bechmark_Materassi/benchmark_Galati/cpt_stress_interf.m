function [sigma_n, sigma_t] = cpt_stress_interf(stress, nodePairsData)
    nni = numel(nodePairsData);
    sigma_n = zeros(nni,1);
    sigma_t = zeros(nni,2);
    for i = 1:nni
        loc_stress = stress(i,:)';
        n = nodePairsData(i).normal;
        t1 = nodePairsData(i).t1;
        t2 = nodePairsData(i).t2;
        N = cpt_normal(n);
        sigma_n(i) = n'*N*loc_stress;
        sigma_t(i,1) = t1'*N*loc_stress;
        sigma_t(i,2) = t2'*N*loc_stress;
    end

end    

















% 
%     nni = numel(nodePairsData);
%     sigma_n = zeros(nni,1);
%     sigma_t = zeros(nni,2);
%     for i = 1:nni
%         ntop = nodePairsData(i).ntop;
%         loc_stress = stress(ntop,:)';
%         n = nodePairsData(i).normal;
%         N = cpt_normal(n);
%         
%         % TODO: check orientation 
%         % Compute two orthonormal tangential basis vectors
%         if abs(n(2)) < 1e-6 && abs(n(3)) < 1e-6  % If n is aligned with x-axis
%             e1 = [0;1;0];  % Pick a different arbitrary vector
%         else
%             e1 = [1;0;0];  
%         end
%         t1 = cross(e1, n);
%         t1 = t1 / norm(t1);
%         t2 = cross(n, t1);
%         t2 = t2 / norm(t2);
% 
%         sigma_n(i) = n'*N*loc_stress;
%         sigma_t(i,1) = t1'*N*loc_stress;
%         sigma_t(i,2) = t2'*N*loc_stress;
%     end
% 
% end
