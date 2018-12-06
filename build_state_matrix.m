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

    for evi = 1:(Nevents_current_ti-1) % loop through the events that occured throughout this track
        % the last frame (starting on the next event) will be over-written in all cases but the last.
        state{ti}( :, SoE(evi,1 ): SoE(evi+1,1) ) = ones( Nsubtracks(ti), SoE(evi+1,1)- SoE(evi,1 ) + 1 ).*(SoE_statemat{ti}(evi,:)');
    end

    % Now over-write any temporary gap-frames in the tracjs with NaN state
    % these are frames where a particle briefly disappears and then reappears
    % these frames are eliminated from consideration.
    state{ti}( isnan(trackmat_xyl(ti).Lamp) ) = NaN ;

end % --- finished for-loop over ti through all non-ephemeral compound track sets.

end % end of function
