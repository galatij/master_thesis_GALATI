function [sigma_n_gp, sigma_t_gp, u_gp] = cpt_gp_quantities(stress_gp, u_gp, interfData, ngauss)
    
    
    ni = numel(interfData);
    sigma_n_gp = cell(ni,1);
    sigma_t_gp = cell(ni,1);
    u_n_gp = 

    for i = 1:ni
        loc_stress_gp = stress_gp{i}; % 4x6 [gauss nodes, components]
        loc_u_gp = u_gp{i};
        n = interfData(i).normal;
        t1 = interfData(i).t1;
        t2 = interfData(i).t2;
        N = cpt_normal(n);
        sigma_n_gp{i} = n'*N*loc_stress_gp';
        sigma_t_gp{i} = t1'*N*loc_stress_gp';
        sigma_t_gp{i} = t2'*N*loc_stress_gp';
        
    end


end    





