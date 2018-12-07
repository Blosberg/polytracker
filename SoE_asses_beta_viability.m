function [ beta_viable ] = SoE_asses_beta_viability( SoE_in, alpha_track, SoE_statemat )
% consider an alpha track, and assess whether its parent could spare an extra particle.

% SeqofEvents (SoE), matrix rows follow this convention:
% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

alpha_birth = find_birth_event(  SoE_in, alpha_track );

beta_track = SoE_in(alpha_birth, 4);

Track_under_focus=beta_track;      % the track we are currently looking at to see if it is carrying away superfluous particles.
event_focus_commence=alpha_birth;  % the point in time from which we begin to look at this particle.
finished = false;

while( ~finished)

   death = find_death_event( SoE_in, Track_under_focus);

   if ( ~all( SoE_statemat( event_focus_commence:(death-1), Track_under_focus) >=2 ) )
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
      event_focus_commence = death;
      Track_under_focus = SoE_in(death,4);
   end

end


end
