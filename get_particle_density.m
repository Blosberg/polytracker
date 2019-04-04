function [  density_mat ] =  get_particle_density( state_matrices_allti, max_state, Nframes, Area )
% produce plots of total track numbers as a function of time



%  monomer_num, track_num, weighted_polymer_num
%  monomer_num           = zeros(1, Nframes);
%  track_num             = zeros(1, Nframes);
%  weighted_polymer_num  = zeros(1, Nframes);


 Num_comp_tracks       = length( state_matrices_allti );
 density_mat           = zeros(max_state, Nframes);


 % Go through all compound track sets.
  for tr = 1: Num_comp_tracks

    % convert NaNs to 0s for this operation
    % @@@ NEED TO SEE WHAT THESE STRUCTURES LOOK LIKE
    state_nansremoved = state_matrices_allti{tr};
    state_nansremoved( isnan(state_nansremoved) )= 0;

    for s = 1:max_state
        density_mat(s,:) = density_mat(s,:) + sum( state_nansremoved == s ,1 );
    end

  end

end
