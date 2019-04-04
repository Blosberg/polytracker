function [ state_lifetime_Plist ] =  get_state_lifetimes ( tracks_input, state_matrices_allti, max_state, dt, Nframes )
% Calculations the number of intervening frames within dimerization events.
% This function grabs lifetimes of each polymer state (i.e. [1] monomer,
% [2] dimer, etc. and outputs a list for each case. (higher-lever functions
% can then concatenate this to a single list irrespective of polymer
% states.
% Lifetime_list{1} is the list of lifetime observations for polymers in the 1 state
% Lifetime_list{2} "" "" in the 2 state, etc.


% --- initialize Plist
Num_comp_tracks     = length( tracks_input );
for s= 1:max_state
    state_lifetime_Plist{s} = [];
end

% loop through all compound track sets.
for ti = 1: Num_comp_tracks

  Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);

  if( size(state_matrices_allti{ti},2) ~= Nframes )
    disp("state matrix size does not correspond to frame number")
    return
  end

  % now loop through sub-tracks within this compound track
  for  sti = 1:Nsubtracks(ti)

     % collect split and merger events in which this track was the enduring party
     Events_this_track = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,4) == sti  ,: );

     event_obit = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,3) == sti  ,: );
     if( size( event_obit,1) ~= 2 )
         disp("not getting birth-death pair")
         return
     end

     birth_track_sti_via_events = event_obit ( event_obit(:,2)==1, :); birth_track_sti_via_events = birth_track_sti_via_events(1);
     death_track_sti_via_events = event_obit ( event_obit(:,2)==2, :); death_track_sti_via_events = death_track_sti_via_events(1);

     birth_track_sti_via_states = min( find(~isnan(state_matrices_allti{ti}(sti,:) )) );
     death_track_sti_via_states = max( find(~isnan(state_matrices_allti{ti}(sti,:) )) );
     % --- special case: ephemeral sub-tracks (i.e. state is NEVER non-NaN).
     % --- set birth AND death equal to birth defined from event.
     if( size( find(~isnan(state_matrices_allti{ti}(sti,:) )), 2) ==0 )
         birth_track_sti_via_states = birth_track_sti_via_events;
         death_track_sti_via_states = birth_track_sti_via_events;
     end

     if( death_track_sti_via_states < tracks_input(ti).seqOfEvents( end, 1 ) )
         death_track_sti_via_states = death_track_sti_via_states +1;
     end

     % N.B.: Whether a track "dies" on the first frame where it begins
     % to be NaN permanently, *OR* on its last non-NaN frame is
     % inconsistent.
     % This cumbersome convention derives from inconsistency in the 
     % original tracksFinal data structure 
     % (which is why sanity checks allow for diff <= 1 ).

     if( birth_track_sti_via_states ~= birth_track_sti_via_events || abs( death_track_sti_via_states - death_track_sti_via_events) > 1 ) 
        disp("FATAL ERROR: Inconsistent birth/death time point between state and event calculation.")
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
                [ state_this_window, excit_duration ]  = get_state( state_matrices_allti{ti}(sti,:), Events_this_track(index,:), Events_this_track(index+1,:) );
                %                                                           ^merger event^     ,     ^split event^
                state_lifetime_Plist{state_this_window} = [ state_lifetime_Plist{state_this_window}, dt * excit_duration ] ;
            end % done cycling through indices with dimerization profile
        end % done "if" checking for the right shape
     end % done "if" checking for more than 1 event

  end % done looping "sti" over Nsubtracks


end % done looping ti over Ncompoundtracks
