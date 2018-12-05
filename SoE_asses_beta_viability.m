function [ beta_viable ] = SoE_asses_beta_viability( SoE_in, alpha_track, SoE_statemat )
% consider an alpha track, and assess whether its parent could spare an extra particle.

% SeqofEvents (SoE), matrix rows follow this convention:
% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

alpha_birth = find_birth_event(alpha_tracks(tr) );

beta_track = SoE_in(alpha_birth, 4);


Track_following=beta_track;
start_following=alpha_birth;
finished = false;

while( ~finished)

   death = find_death_event(Track_following);

   if ( ~all( SoE_statemat( start_following:(death-1), Track_following) >=2 ) )
      % cannot spare an extra particle before death.
      beta_viable = false;
      finished    = true;
      break;
   elseif( isnan( SoE_in(death,4)) )
      % CAN spare, and dies into nothing.
      beta_viable = true;
      finished    = true;
      break;
   else
      % CAN spare, but dies merging into something else. Follow THAT one
      % now.
      start_following = death;
      Track_following = SoE_in(death,4);
   end

end


end
