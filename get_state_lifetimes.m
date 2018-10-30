function [ state_lifetime_list ] =  get_state_lifetimes ( state_matrix_allti, tracks_input, max_state  )
% Added on boztower at 21:50 on 25.10.2018
% 

Num_comp_tracks = length( tracks_input )


% go through all non-ephemeral compound track sets.
for ti = 1: Num_comp_tracks 
   
  Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);
  
  for  sti = 1:Nsubtracks(ti)
     
     % collect split and merger events in which this track was the enduring party     
     Events_this_track = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,4) == sti  ,: ) 

     if( dim( Events_this_track,1) >1 )
     % look for death events immediately followed by birth events

        indices_right_shape = find (diff( Events_this_track(:,2) )  >= 1 ) 
        % since diff returns an array shorter by one, we can always
        % access the indices in this array, AND the one after.

        for j=1:length( indices_right_shape )

            [state, excit_duration ]   = get_state( state_matrix_allti{ti}(sti,:), Events_this_track(j,:), Events_this_track(j+1,:) )
            state_lifetime_list{state} = [ state_lifetime_list{state}, excit_duration ] 
        end     

     end  
   
  end


end
