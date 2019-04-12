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
    return;
end
% -----

for j= 1:length(indices_right_shape)
    omit_event = false;
    index = indices_right_shape(j);

    % DOUBLE COUNTING CHECK:

    if( Event_list( index, 2 ) == 1 )
        % i.e. if the first of the two events is a "birth/split" (meaning
        % this is an "off-window" we are looking at.
        if ( all( Event_list(index, 3:4) ==  Event_list(index+1, 3:4) ) || all( Event_list(index, 3:4) ==  fliplr( Event_list(index+1, 3:4) ) )  )
            % i.e. For track X, born out of Y, either:
            % X dies by merging immediately back into Y or
            % Y dies by merging immediately back into X
            % with no intermediate interactions.
            % Otherwise, skip this step:

            if ( this_subtr ~= Event_list(index+1, 4) )
                % If the above condition _is_ satisfied, then we only count
                % this "off-window" once --when "this_subtr[ack] is the one
                % that survives the re-merging event.
                % Otherwise, omit this off-windoe to avoid double-counting.
                omit_event = true;
            end
        end

    end

    if ( ~omit_event )

        window_duration = Event_list(index+1, 1) - Event_list(index, 1) + 1;
        onoff_windows   = [ onoff_windows ,  (dt * window_duration )  ] ;

    end

end % done cycling through indices with "on"-shape

end
