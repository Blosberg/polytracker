function [ state, max_state ] =  build_state_matrix ( tracks_input, trackmat_xyl, Nframes )
% build a matrix for each composite track where the integer elements characterize the state of the polymer in that frame.
% trackdat_xyl (built from previous function) is used to sort out missed frames.

max_state = 1;
Num_comp_tracks = length( tracks_input);

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.

    Nevents_current_ti   = size(tracks_input(ti).seqOfEvents,      1);
    Nsubtracks(ti)       = size(tracks_input(ti).tracksCoordAmpCG, 1);

    % initial "adjustment matrix" is all zeros.
    adj_mat_init         = zeros(Nevents_current_ti, Nsubtracks(ti));

    % and the initial state_matrix is calculated according to default assumptions.
    SoE_statemat{ti}     = SoE_time_evolve_statemat ( tracks_input(ti).seqOfEvents, Nsubtracks(ti), adj_mat_init);

    % if states are untenable, then these assumptions are revised.
    if( ~check_states_tenable ( SoE_statemat{ti} ) )
        [ SoE_statemat{ti}, adj_mat{ti} ]= SoE_revise_statemat ( tracks_input(ti).seqOfEvents, SoE_statemat{ti} );
    end

    max_state = max( max_state, max( max( SoE_statemat{ti})) );

end

% =============================================================================
% --- POPULATE THE PER-FRAME STATE MATRIX ----

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.

    SoE                  = tracks_input(ti).seqOfEvents;

    Nsubtracks(ti)       = size( tracks_input(ti).tracksCoordAmpCG, 1);
    Nevents_current_ti   = size( SoE, 1 );

    state{ti}      = NaN( [ Nsubtracks(ti), Nframes ] );

    for evi = 1:(Nevents_current_ti-1) 
    % loop through the events that occured throughout this track
        current_event_frame = SoE(evi  , 1 );
        next_event_frame    = SoE(evi+1, 1 );
        
        
        if( (next_event_frame - current_event_frame) >= 1 )
            % Define state for the current frame up to BUT NOT INCLUDING the 
            % frame of the next event.
            state{ti}( :, current_event_frame : next_event_frame -1 ) = ones( Nsubtracks(ti), next_event_frame - current_event_frame ).*(SoE_statemat{ti}(evi,:)');
        end
    
        simultaneous_events =  find( SoE(:,1)==current_event_frame );
        % over-write states for this FRAME among active tracks (including
        % the present one).
        for sevi = 1:length(simultaneous_events)
            sev = simultaneous_events( sevi ); % sev is the simultaneous event index in the original SoE 
            state{ti}( SoE(sev, 3), SoE(sev, 1) ) = SoE_get_active_track_state( trackmat_xyl(ti).Lamp( SoE(sev,3), SoE(sev,1) ), SoE_statemat{ti}, sev, SoE(sev,3) );
        end
        %                                                                                          ^ active track, ^ time                               
    end

    % now the last one:
    evi = Nevents_current_ti;
    current_event_frame = SoE(Nevents_current_ti, 1 );    

    % set all non-active tracks.
    state{ti}( :, SoE(evi,1 ) )          = ones( Nsubtracks(ti), 1 ).*(SoE_statemat{ti}(evi,:)');

    simultaneous_events =  find( SoE(:,1)==current_event_frame );
    % over-write states for this FRAME among active tracks (including
    % the present one).
    for sevi = 1:length(simultaneous_events)
        sev = simultaneous_events( sevi ); % sev is the simultaneous event index in the original SoE
        state{ti}( SoE(sev, 3), SoE(sev, 1) ) = SoE_get_active_track_state( trackmat_xyl(ti).Lamp( SoE(sev,3), SoE(sev,1) ), SoE_statemat{ti}, sev, SoE(sev,3) );
    end
    
    % We're now through to the end of the movie for this compound track.
    % Last step: over-write any gap-frames (i.e. transient frames where the
    % particle disappeared, and then reappeard after). Set these states to
    % NaN.    
    state{ti}( isnan(trackmat_xyl(ti).Lamp)  ) = NaN;
    
end % --- finished for-loop over ti through all non-ephemeral compound track sets.

end % end of function
