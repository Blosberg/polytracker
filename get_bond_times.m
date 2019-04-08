function [ bond_times_on, bond_times_off ] =  get_bond_times ( tracks_input, state_matrices_allti, max_state, dt, Nframes )
% Calculations the number of intervening frames within dimerization events.
% unlike previous script, this is state-agnostic. Just on/off times
% regardless of state.


bond_times_on  = [];
bond_times_off = [];

Ncomptracks  = size(tracks_input);

% loop through all compound track sets.
for tr = 1: Ncomptracks

  Nsubtracks(tr) =  size( tracks_input(tr).tracksCoordAmpCG, 1);

  if( size(state_matrices_allti{tr},2) ~= Nframes )
    disp("state matrix size does not correspond to frame number")
    return
  end

  % now loop through sub-tracks within this compound track
  for  subtr = 1:Nsubtracks(tr)

     % collect split and merger events in which this track (subtr) was the enduring party
     Events_this_track_enduring = tracks_input(tr).seqOfEvents( tracks_input(tr).seqOfEvents(:,4) == subtr  ,: );

     % collect split and merger events in which this track (subtr) was somehow involved at all
     Events_this_track_involved = tracks_input(tr).seqOfEvents( (tracks_input(tr).seqOfEvents(:,4) == subtr | tracks_input(tr).seqOfEvents(:,3) == subtr ),: );

     % Filter only events for which the secondary interacting track is not "NaN".
     % i.e. we want split/merger events, not birth/death out of nothing.
     Events_this_track_involved = Events_this_track_involved( ~isnan(Events_this_track_involved(:,4)), : );

     % collect _just_ the birth and death events (the obit[uary]) of this
     % subtrack (subtr) to determine the lifetime.
     event_obit = tracks_input(tr).seqOfEvents( tracks_input(tr).seqOfEvents(:,3) == subtr  ,: );
     if( size( event_obit,1) ~= 2 )
         disp("not getting birth-death pair")
         return
     end

     % @@@ TODO: Should just be able to take first and second elements, no?
     birth_track_sti_via_events = event_obit ( event_obit(:,2)==1, :); birth_track_sti_via_events = birth_track_sti_via_events(1);
     death_track_sti_via_events = event_obit ( event_obit(:,2)==2, :); death_track_sti_via_events = death_track_sti_via_events(1);

     birth_track_sti_via_states = min( find(~isnan(state_matrices_allti{tr}(subtr,:) )) );
     death_track_sti_via_states = max( find(~isnan(state_matrices_allti{tr}(subtr,:) )) );
     % --- special case: ephemeral sub-tracks (i.e. state is NEVER non-NaN).
     % --- set birth AND death equal to birth defined from event.
     if( size( find(~isnan(state_matrices_allti{tr}(subtr,:) )), 2) ==0 )
         birth_track_sti_via_states = birth_track_sti_via_events;
         death_track_sti_via_states = birth_track_sti_via_events;
     end

     if( death_track_sti_via_states < tracks_input(tr).seqOfEvents( end, 1 ) )
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



   % an "on" window is when this track endures a merge followed directly by
   % a split ( birth - death = -1 )
   bond_times_on  = [bond_times_on  , get_onoff_windows(Events_this_track_enduring, -1, dt, subtr) ];

   % An "off" window is when this track is one of the tracks involved in a
   % split, and immediately thereafter is one of the tracks involved in a
   % merge: (death - birth = 1 )
   bond_times_off = [bond_times_off , get_onoff_windows(Events_this_track_involved, 1,  dt, subtr) ];



  end % done looping "sti" over Nsubtracks


end % done looping ti over Ncompoundtracks
