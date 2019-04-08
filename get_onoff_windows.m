function [ onoff_windows ] =  get_onoff_windows( Event_list, bddiff_target, dt, this_subtr )
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
% -----

% If none are the right "shape" then we are finished.
if( length(indices_right_shape) < 1)
    return
end
% -----

for j= 1:length(indices_right_shape)
    omit_event = false;
    index = indices_right_shape(j);

    % If the same two particles split and then merge again, we have to avoid
    % double-counting this event:
    if( Event_list( index, 2 ) == 1 )
        if ( Event_list(index, 3:4) ==  Event_list(index+1, 3:4) )
            disp("Case 1: track dies from whence it was born")
            if ( this_subtr ~= Event_list(index+1, 4) )
                omit_event = true;
            end
        else ( Event_list(index, 3:4) ==  fliplr( Event_list(index+1, 3:4) ) )
            disp("Case 2: parent dies into the child")
            if ( this_subtr ~= Event_list(index+1, 4) )
                omit_event = true;
            end
        end

    end

    if ( ~omit_event )

        window_on = Event_list(index+1, 1) - Event_list(index, 1) + 1;
        onoff_windows = [ onoff_windows ,  (dt * window_on )  ] ;

    end

end % done cycling through indices with "on"-shape

end
