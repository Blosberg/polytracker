function [ state, max_state ] =  build_state_matrix ( tracks_input, trackmat_xyl, Nframes )
% build a matrix for each composite track where the integer elements characterize the state of the polymer in that frame.
% trackdat_xyl (built from previous function) is used to sort out missed frames.

max_state = 1;
Num_comp_tracks = length( tracks_input);

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.

    Nevents    = size(tracks_input(ti).seqOfEvents,1);
    Nsubtracks = size(tracks_input(ti).tracksCoordAmpCG,1);

    % initial "adjustment matrix" is all zeros.
    adj_mat_init     = zeros(Nevents, Nsubtracks);

    % and the initial state_matrix is calculated according to default assumptions.
    SoE_statemat{ti} = SoE_time_evolve_matrix ( tracks_input(ti).seqOfEvents,Nsubtracks, adj_mat_init);

    % if states are untenable, then these assumptions are revised.
    if( ~check_states_tenable ( SoE_statemat{ti} ) )
        [ SoE_statemat{ti}, adj_mat{ti} ]= SoE_revise_statmat ( tracks_input(ti).seqOfEvents, SoE_statemat{ti} )
    end
end

% =============================================================================
% --- POPULATE THE PER-FRAME STATE MATRIX ----

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.

    SoE            = tracks_input(ti).seqOfEvents;

    Nsubtracks(ti) = size( tracks_input(ti).tracksCoordAmpCG, 1);
    Nevents        = size( SoE );

    state{ti}      = NaNs( Nsubtracks(ti), Nframes );

    for evi = 1:(Nevents-1) % loop through the events that occured throughout this track

       if( SoE(evi+1,1)- SoE(evi,1 ) > 0 )
           % check that this event lasts at least one frame:
           state{ti}( :, SoE(evi,1 ):SoE(evi+1,1)-1 ) = ones( Nsubtracks, SoE(evi+1,1)- SoE(evi,1 ) ).*SoE_statemat{ti}(evi,:)';
       end

    end
    state{ti}( :, SoE(evi,1 ):end) = ones( Nsubtracks, Nframes - SoE(evi,1 ) ).*SoE_statemat{ti}(Nevents,:)';

end


end
