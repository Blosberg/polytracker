function [ SoE_state_matrix ] =  SoE_time_evolve_statemat ( SoE, Nsubtracks, adj_mat )
% Input: matrix sequence of events
% Output: matrix of tracks vs. Events. Uncorrected. values may be 0 or negative,
% and, if so, will be corrected by a separate script.

Nevents               = size( SoE, 1 );
SoE_state_matrix_init = NaN( [ Nevents, Nsubtracks] );

% initialize:
prev_states  = NaN( [ 1, Nsubtracks ] );

% cycle through events:
for evi = 1:Nevents

    if( evi > 1 )
       prev_states = SoE_state_matrix_init(evi-1,:);
    end

    % preserve states from previous event for tracks irrelevant to the present event
    % ( relevant tracks will be over-written later)
    SoE_state_matrix_init( evi, : )          = prev_states;

    %-------------------- BIRTH: --------------------
    if ( SoE(evi,2) == 1)
       SoE_state_matrix_init( evi, SoE(evi,3) ) = 1 + adj_mat(evi, SoE(evi,3) ) ; 
       %   index for this track:  ^ 
       % initialize as monomer + whatever is indicated by the adjustment matrix

       % ---- born through split?
       if ( ~isnan( SoE(evi,4) )  )
           SoE_state_matrix_init( evi, SoE(evi,4) ) = prev_states ( SoE(evi, 4) ) -1 - adj_mat(evi, SoE(evi,3)) ; 
       %   index for other track:     ^                             ^ 
       % Parent track decreases by (1+adjustment_matrix of born track at this event), 
       % The default assumption being that the born track is always is a 
       %  monomer state). 
       end

    %-------------------- DEATH: --------------------
    else if ( SoE(evi,2) == 2 )

       SoE_state_matrix_init( evi, SoE(evi,3) ) = NaN;  
       %   index for this track:  ^ 

       % ---- death through merger ?
       if ( ~isnan( SoE(evi,4) )  )
           SoE_state_matrix_init( evi, SoE(evi,4) ) = prev_states ( SoE(evi, 4) ) + prev_states ( SoE(evi, 3) ) ; % increment by state of dying track
       %   index for other track:     ^ 
       end

    end

end % end of loop through Events.


end % end of function
