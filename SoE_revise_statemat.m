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

    % "alpha" = tracks that merge into broken track between the latter's birth and breakage, AND which were born from split events
    % (from which they could potentially bring an additional particle)
    [ alpha_tracks, alpha_events ] = SoE_get_alpha_merge_list( SoE, broken_track, breaking_event);


    % initialize m to the last merge event
    % and then work backwards:

    m = length( alpha_tracks);
    while( m > 0  && ~adj_mat_changed )

       beta_viable = SoE_asses_beta_viability( SoE, alpha_tracks(m), SoE_statemat_out)

       if(beta_viable)
          alpha_birth = find_birth_event(alpha_tracks(m));
          adj_mat( alpha_birth, alpha_tracks(m) ) = adj_mat( alpha_birth, alpha_tracks(m) ) +1
          adj_mat_changed = true;
          break;
       else
          m = m-1;
       end
    end

    % if we still haven't found a
    % candidate: then increment the broken track at the outset.
    if ( ~adj_mat_changed )
        adj_mat( birthevent_broken, broken_track ) = adj_mat( birthevent_broken, broken_track ) + 1;
        adj_mat_changed = true;
    end

    % rebuild the SoE_statemat with the new adj_mat:
    SoE_statemat_out = time_evolve_SoE_statemat ( SoE, Nsubtracks, adj_mat )
end % end of loop "while states are untenable".

end % end of function
