function [ alpha_merge_list ] = SoE_get_alpha_merge_list( SoE_in, Track, Event )
% get the list of tracks that died through merging into "Track" before "Event"
% SeqofEvents (SoE), matrix rows follow this convention:
% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

alpha_merge_list = [];

if ( Event == 1 )
   return; % @@@ TODO: check that the right value is returned here
else
   SoE_timely = SoE_in(1:(Event-1), : );
end

% list of events for deaths by merging into "Track"
death_merger_events = find( (SoE_timely(:,2) == 2).*(SoE_timely(:,4)==Track) );

if ( length(death_merger_events) >= 1 )
    
    % corresponding track numbers:
    merging_track_list = SoE_in( death_merger_events, 3);
    
    for tr = 1:length( merging_track_list )
      
       alpha_birth_event = find_birth_event(SoE_in, merging_track_list(tr));
       if( ~isnan(SoE_in(alpha_birth_event,4)) )
          alpha_merge_list = [alpha_merge_list, merging_track_list(tr)];
       end

    end
end

% alpha_merge_list is now the list of all tracks that die through merging
% with "Track"  before the latter becomes untenable, and are born through
% splitting off of another track. Most of the time, this will simply be
% empty.

end
