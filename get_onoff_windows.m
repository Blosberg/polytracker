function [ onoff_windows ] =  get_onoff_windows( Event_list, bddiff_target, dt )
% identify windows of time in which the dimerization is either "on" or
% "off" for building histograms.
% Event_list is a list of events that involve the track of interest in the
% right way (depending on whether we want an "on" window or an "off" window
% bddiff_target is either 1 or -1 if we want on or off window respetively.

onoff_windows = [];

% only consider this tracklist if there are at least 2 events.
if( size( Event_list,1) <= 1 )
    return;
end



%  look for birth/death diff with the right features.
indices_right_shape = find (diff( Event_list(:,2) )  == bddiff_target );

% since diff returns an array shorter by one, we can always
% access the indices in this array, AND the one after.

% If none are the right "shape" then we are finished.
if( length(indices_right_shape) < 1)
    return
end


for j= 1:length(indices_right_shape)
    index = indices_right_shape(j);

    window_on = Event_list(index+1, 1) - Event_list(index, 1) + 1;

    onoff_windows = [ onoff_windows ,  (dt * window_on )  ] ;
end % done cycling through indices with "on"-shape



end
