function [ state_lifetime_list ] =  get_state_lifetimes ( tracks_input, state_matrix_allti, max_state , Nframes )
% Added on boztower at 21:50 on 25.10.2018
%

Num_comp_tracks     = length( tracks_input );
for s= 1:max_state
    state_lifetime_list{s} = [];
end

% go through all non-ephemeral compound track sets.
for ti = 1: Num_comp_tracks

  % disp(ti)
  Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);

  if( size(state_matrix_allti{ti},2) ~= Nframes )
    disp("state matrix size does not correspond to frame number")
    return
  end


  for  sti = 1:Nsubtracks(ti)

     % collect split and merger events in which this track was the enduring party
     Events_this_track = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,4) == sti  ,: );

     state_birth = min( find(~isnan(state_matrix_allti{ti}(sti,:) )) );
     state_death = max( find(~isnan(state_matrix_allti{ti}(sti,:) )) ) +1;
     % track "dies" on the first frame where it begins to be NaN permanently.
     % even if this is a "hypothetical" frame, immediately following the end of the movie.

     % events that list the birth and death of this subtrack
     event_obit = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,3) == sti  ,: );
     if( size( event_obit,1) ~= 2 )
         disp("not getting birth-death pair")
         return
     end

     event_birth = event_obit ( event_obit(:,2)==1, :); event_birth = event_birth(1);
     event_death = event_obit ( event_obit(:,2)==2, :); event_death = event_death(1);

     % terminal event? if so, then it "dies" on what _would have_ been the next frame
     % original developpers were EXTREMELY sloppy on this point
     % By their notation, the frame # from seqOfEvents points to
     % ---- the last frame before event, if either (A) the event was terminal to the compound track
     %                                          OR (B) Was the same frame as the birth.
     %                                          OR (C) some other as-yet undetermined condition (!!)
     %  if none of these things were true, then USUALLY, the event
     % refers to first post-event frame.
     % -- A convoluted and imprecise convention.
     %  I've tried my best to ensure that state-changes are ALWAYS noted immediately
     %  after their causes.

     if ( event_death  == max( tracks_input(ti).seqOfEvents(:,1) ) || event_death == event_birth)
       event_death = event_death +1;
     end

     % and sometimes the event is recorded before-frame anyway, for no apparent reason.
     if( state_birth ~= event_birth || abs( state_death - event_death)>1 )
        disp("Inconsistent birth/death time point between state and event calculation.")
        return
     end

     if( size( Events_this_track,1) >1 )
     % look for death events immediately followed by birth events

        indices_right_shape = find (diff( Events_this_track(:,2) )  <= -1 );
        % since diff returns an array shorter by one, we can always
        % access the indices in this array, AND the one after.

        if( length(indices_right_shape) >= 1)

            for j= 1:length(indices_right_shape)
                index = indices_right_shape(j);
                [ state_this_window, excit_duration ]  = get_state( state_matrix_allti{ti}(sti,:), Events_this_track(index,:), Events_this_track(index+1,:) );
                %                                                           ^merger event^     ,     ^split event^
                state_lifetime_list{state_this_window} = [ state_lifetime_list{state_this_window}, excit_duration ] ;
            end
        end
     end

  end


end

% ===================================================
% Following documentation taken from the comments of plotComptrack:
% Within each Tracks(i) data, there will be the following elements:

%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.


