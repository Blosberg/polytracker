function [ active_track_state ] = SoE_get_active_track_state( Lamp_current_frame, SoE_statemat, evi, tracknum )
%Gathers a state vector for a given frame

if ( isnan(Lamp_current_frame) )
    % active track is currently illuminating "NaN" ==> it must be dead.
    active_track_state = NaN;
else
    % Track is still actively illuminating for this frame. It must be live.
    % therefore, its state must be a finite value.

    % if we have a finite state value for the onset of this frame, use it.
    if( ~isnan( SoE_statemat(evi, tracknum))  )
        active_track_state =  SoE_statemat(evi, tracknum);
    else
        % If state just became NaN, look at the previous state.
        if( isnan( SoE_statemat(evi-1, tracknum)  ) )
            % if that is ALSO NaN, then flag this error and halt.
            disp("ERROR: cannot find any valid state.");
            return;
        else
            % otherwise, assume that this Event is the terminating frame,
            % and preserve the last state recorded.
            active_track_state = SoE_statemat(evi-1, tracknum);
        end        
        
    end % end of "Is there a finite state for the current event."
    
end % end of "is the track still Illuminating".

end % end of function
