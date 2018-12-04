function [ SoE_statemat_out, adj_mat ] =  SoE_revise_statemat ( SoE, SoE_statemat_in )
% Input: matrix sequence of events, matrix of states

% SeqofEvents (SoE), matrix rows follow this convention:
% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

Nevents      = size( SoE, 1 );
Nsubtracks   = size( SoE_statemat, 2 );

% initial "adjustment matrix" is all zeros.
adj_mat      = zeros(Nevents, Nsubtracks);


% initialize:
SoE_statemat_out  = SoE_statemat_in;

% continue through this loop until check_states_tenable returns 1. Then we have our result
while ( ~check_states_tenable (SoE_statemat_out)   )

    % has the adjustment matrix changed in this iteration?
    adj_mat_changed = false;

    % Identify the track for which initial untenability begins *last*
    % and the Event at which it first becomes untenable
    [ broken_track, breaking_event ] = find_where_broken( SoE_statemat_out )

    % find event of this tracks birth
    birthevent_broken = find_birth_event( SoE, broken_track );

    % "alpha" = tracks that merge into broken track between the latter's birth and breakage
    [ alpha_tracks, alpha_events ] = SoE_get_incoming_merge_list( SoE, broken_track, breaking_event);

    M = 0;

    if( length( alpha_tracks) > 0 )
       beta_tracks = 
    end

    %initialize m to the last merge event
    m = M;

    while( m > 0 ) 

       @@@ if i
       ... adj_mat_changed = true;
       m = m-1;
    end 

    % if we get to the end of this list and still haven't found a suitable candidate:
    % then embiggen (sure that's a word) the broken track at the outset.
    if ( ~adj_mat_changed ) 
        adj_mat( birthevent_broken, broken_track );
        adj_mat_changed = true;
    end

    % rebuild the SoE_statemat with the new adj_mat:
    SoE_statemat_out = time_evolve_SoE_statemat ( SoE, Nsubtracks, adj_mat )
end % end of loop through Events.




end % end of function
