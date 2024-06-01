function rhs = assemble_rhs_neu(nn,time_cur,neum_force);
%-----------------------------------------------------------------------------------------
%
% Uses the neum_force structure to build the rhs. The force at current time is
% interpolated from the forces at the previous and subsequent ones.
%
%-----------------------------------------------------------------------------------------

times = neum_force.times;
nstep = numel(times);
ind_prev = find(times>time_cur,1,'first')-1;
ind_next = ind_prev + 1;
if ind_prev < 0 || ind_next > nstep
   error('Time outside possible range');
end
time_prev = times(ind_prev);
time_next = times(ind_next);
fac = (time_cur-time_prev) / (time_next-time_prev);
forces = neum_force.forces;
ff = forces(:,ind_prev) + fac*(forces(:,ind_next) - forces(:,ind_prev));

% Assemble rhs
rhs = zeros(nn,1);
ind = neum_force.ind_neu;
rhs(ind) = ff;

end
