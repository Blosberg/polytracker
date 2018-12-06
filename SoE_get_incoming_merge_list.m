function [ merging_tracks, merge_events ] = SoE_get_incoming_merge_list( SoE, broken_track, breaking_event)
% get the list event tracks and event indices for mergers into the broken track after it is born, but before breaking_event.

% find event of this tracks birth
birth = find( SoE(:,2) == 1 && SoE(:,3) == broken_track );

if ( length( birth) ~= 1 )
   disp("ERROR: non-unique birth event for broken_track");
end

% identify events where another track died by merging into this one
merge_events  = find( SoE(:,2)==2 && SoE(:,4)==broken_track );

% take the subset of those that occured before breakage:
merge_events  = merge_events( merge_events < breaking_event );

if ( length( find( merge_events < birth ) ) >= 1  )
   disp("ERROR: some tracks are merging before this track is born. Investigate!")
end

merging_tracks = SoE( merge_events, 3);

end
