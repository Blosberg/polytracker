function [ SoE_state_matrix_revised ] =  build_SoE_state_matrix_init ( SoE )
% Input: matrix sequence of events


% SeqofEvents (SoE), matrix rows follow this convention:
% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

Nevents      = size( SoE, 1 );
State_Ev_mat = NaN( [ Nevents, Nsubtracks] );


% initialize:
prev_states  = NaN( [ 1, Nsubtracks ] );

% cycle through events:
for evi = 1:Nevents

    if( evi > 1 )
       prev_states = State_ev_mat(evi-1,:);
    end

    % set all states equal to previous event ( relevant tracks will be over-written)
    State_Ev_mat( evi, : )          = prev_states;

    %-------------------- BIRTH: --------------------
    if ( SoE(evi,2) == 1)
       State_Ev_mat( evi, SoE(evi,3) ) = 1; % initialize as monomer
       %   index for this track:  ^ 

       % ---- born through split?
       if ( ~isnan( SoE(evi,4) )  )
           State_Ev_mat( evi, SoE(evi,4) ) = prev_states ( SoE(evi, 4) ) -1; 
       %   index for other track:     ^                             ^ 
       % Parent track is decremented by 1, under the assumption that the born
       % track is always initialized in the monomer state) This assumption is
       % revisited only when necessary in another function. 
       end

    %-------------------- DEATH: --------------------
    else if ( SoE(evi,2) == 2 )

       State_Ev_mat( evi, SoE(evi,3) ) = NaN;  
       %   index for this track:  ^ 

       % ---- death through merger ?
       if ( ~isnan( SoE(evi,4) )  )
           State_Ev_mat( evi, SoE(evi,4) ) = prev_states ( SoE(evi, 4) ) + prev_states ( SoE(evi, 3) ) ; % increment by state of dying track
       %   index for other track:     ^ 
       end

    end

end % end of loop through Events.




end % end of function
