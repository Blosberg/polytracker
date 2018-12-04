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
% 
% =============================================================================
% --- POPULATE THE PER-FRAME STATE MATRIX ----

for ti = 1:Num_comp_tracks % go through all non-ephemeral independent track sets.

    Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);  
      
    state{ti} = ones( Nsubtracks(ti), Nframes );

    %initialize "state" by assuming monomer status throughout entire movie.
   
    for evi = 1:size( tracks_input(ti).seqOfEvents, 1) % loop through the events that occured throughout this track
       
        if ( tracks_input(ti).seqOfEvents(evi,2) == 1)% ------- BIRTH:
           
            % ---- set states for _THIS_ subtrack. Is it delayed from 1 ?
            if ( tracks_input(ti).seqOfEvents(evi,1) > 1  )
                % subtrack does not begin at frame 1; set state for all previous time-points to NaN
                state{ti}( tracks_input(ti).seqOfEvents(evi,3), 1:((tracks_input(ti).seqOfEvents(evi,1))-1)) = NaN ;
                %                    ^ col. 3 refers to  the subtrack that just got born.
            end
           
            % ---- If born through split, update state of corresponding track
            if ( ~isnan( tracks_input(ti).seqOfEvents(evi,4) )  )
                % if so, then decrement state of the other track affected by this
                % point on to the end of the full track.
                temp_dec = state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), tracks_input(ti).seqOfEvents(evi,1) : Nframes );
                           state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), tracks_input(ti).seqOfEvents(evi,1) : Nframes ) = temp_dec - 1;
                %                    |                                      |    |                                  |              |             |
                %                    |^row number corresponding to index of |    | ^ Frame  # provided by 1rst col. |              |             |
                %                    |other subtrack affected by this event |    | From this event to end time      |              | decrement state |
            end
        end
       
        if ( tracks_input(ti).seqOfEvents(evi,2) == 2 )% ------- DEATH:
            % evi describes a Track that has just "died":
           
            % ---- non-terminal case?
            if ( tracks_input(ti).seqOfEvents(evi,1) < Nframes  )
                % subtrack does not end at the very last frame; set all states
                % from this point to the end as NaN
                state{ti}( tracks_input(ti).seqOfEvents(evi,3) , ((tracks_input(ti).seqOfEvents(evi,1))+1): Nframes ) = NaN ;
                %                                           ^ col. 3 =subtrack that just died.
                       
             end % finished "if" this event was non-terminal (determines 
          
                          % ---- death through merger ?
             if ( ~isnan( tracks_input(ti).seqOfEvents(evi,4) )  )
                   % if so, then increment state of the other track affected by this
                   % point on to the end of the full track.
                   temp_inc = state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), (tracks_input(ti).seqOfEvents(evi,1) ): Nframes );
                              state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), (tracks_input(ti).seqOfEvents(evi,1) ): Nframes ) = temp_inc + 1 ;
                   %                    |                                  |   |                                            |    |                 |
                   %                    |^row number == index of other     |   | "now" == Frame (from 1rst col.)^           |    |                 |
                   %                    | subtrack affected by this event  |   |             from now to end of run   ^     |    | increment state |
               
             end % --- finished "if" checking for merger
          
        end % --- finished "if" checking for death
            
       
    end % --- finished "for loop" over evi
  

    %---- set states relative to baseline (minimal == 1) 
    for sti = 1:Nsubtracks(ti) % within each subtrack set minimal finite state to 1, and all others relative to that.
        offset = 1 - min( state{ti}(sti,:) );
        state{ti}(sti,:) = state{ti}(sti,:) + offset; % ==   
        
        state{ti}(sti,  isnan(trackmat_xyl(ti).Lamp(sti,:)) ) = NaN ;
    end
   
    % max_state => largest state seen so far
    max_state = max( [ max( state{ti} ), max_state  ]);
   
end % --- finished for-loop over ti through all non-ephemeral independent track sets.

end
