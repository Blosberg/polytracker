function [ tracknum_density, weighted_polymer_dens ] =  get_particle_density( state, Nframes )
% produce plots of total track numbers as a function of time

 Num_comp_tracks     = length( state );

 tracknum_density       = zeros(1, Nframes); 
 weighted_polymer_dens  = zeros(1, Nframes);

 % go through all non-ephemeral compound track sets.
  for ti = 1: Num_comp_tracks 


 % convert NaNs to 0s for this operation
 state_zerod = state{ti};
 state_zerod( isnan(state_zerod) )= 0;
 
 weighted_polymer_dens  =  weighted_polymer_dens + sum( state_zerod, 1 );
 
 state_oned = sum( (state_zerod > 0), 1);
 
 tracknum_density = tracknum_density + state_oned;
 end 

end
