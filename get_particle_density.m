function [  monomer_num, track_num, weighted_polymer_num ] =  get_particle_density( state, Nframes )
% produce plots of total track numbers as a function of time

 Num_comp_tracks       = length( state );

 monomer_num           = zeros(1, Nframes);
 track_num             = zeros(1, Nframes); 
 weighted_polymer_num  = zeros(1, Nframes);

 % Go through all compound track sets.
  for ti = 1: Num_comp_tracks 

    % convert NaNs to 0s for this operation
    state_nansremoved = state{ti};
    state_nansremoved( isnan(state_nansremoved) )= 0;

    % cumulative number of particles:
    weighted_polymer_num  =  weighted_polymer_num + sum( state_nansremoved, 1 );
 
    % number of monomers is sum of wherever states = 1
    state_monod   = sum( (state_nansremoved == 1), 1);
    monomer_num   = monomer_num  + state_monod;

    % number of tracks (regardless of state): 
    state_notzero = sum( (state_nansremoved > 0), 1);
    track_num = track_num + state_notzero;
    
  end

end
