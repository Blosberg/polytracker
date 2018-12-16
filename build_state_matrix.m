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

    evi = 1;
    prev_states = SoE_state_matrix(1,:);

    while( evi < Nevents_current_ti)
    % loop through the events that occured throughout this track
        current_event_frame = SoE(evi  , 1 );

        % if multiple events happen simultaneously, use the states after the last of them.
        simultaneous_events = find( SoE(:,1)==current_event_frame );
        evi                 = max(simultaneous_events);

        if( evi  == size(SoE,1)
            % if the movie is ending here, then assign the states of this one frame:
            state{ti}( :, current_event_frame ) = ones( Nsubtracks(ti), 1 ).*(prev_states');
        else
	    % otherwise Define state for the current frame up to BUT NOT
	    % INCLUDING the frame of the next event.

            next_event_frame    = SoE(evi+1, 1 );
            state{ti}( :, current_event_frame : next_event_frame -1 ) = ones( Nsubtracks(ti), next_event_frame - current_event_frame ).*(SoE_statemat{ti}(evi,:)');

            % set "prev_states" for next iteration
            prev_states = SoE_statemat{ti}(evi,:);
        end

        evi = evi + 1;

    end % terminate while loop over event index "evi"

    % We're now through to the end of the movie for this compound track.
    % Last step: over-write any gap-frames (i.e. transient frames where the
    % particle disappeared, and then reappeard after). Set these states to
    % NaN.
    state{ti}( isnan(trackmat_xyl(ti).Lamp)  ) = NaN;

end % --- finished for-loop over ti through all non-ephemeral compound track sets.

% ======== DOUBLE CHECK : ===============

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.
  for  sti = 1:Nsubtracks(ti)

     SoE                        = tracks_input(ti).seqOfEvents;
     event_obit                 = SoE( SoE(:,3) == sti , : );
     birth_track_sti_via_events = event_obit ( event_obit(:,2)==1, :); birth_track_sti_via_events = birth_track_sti_via_events(1);
     death_track_sti_via_events = event_obit ( event_obit(:,2)==2, :); death_track_sti_via_events = death_track_sti_via_events(1);


     birth_track_sti_via_states = min( find(~isnan(state{ti}(sti,:) )) );
     death_track_sti_via_states = max( find(~isnan(state{ti}(sti,:) )) );
     if( death_track_sti_via_states < tracks_input(ti).seqOfEvents( end, 1 );
         death_track_sti_via_states = death_track_sti_via_states +1;
     end

     if( birth_track_sti_via_states ~= birth_track_sti_via_events || abs( death_track_sti_via_states - death_track_sti_via_events) > 0 ) % @@@ may need to change this back to >1
        disp("Inconsistent birth/death time point between state and event calculation.")
        return
     end



  end
end

% ====== FINISHED DOUBLE CHECK : =======


end % end of function
