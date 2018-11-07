  function [ tracksdat, Nframes ] =  purge_ephemeral ( tracks_input_RAW )
    
  
  Num_indep_tracks_RAW = length(tracks_input_RAW);
  Nframes = 1;
  
  for ti = 1:Num_indep_tracks_RAW
      
      % Nframes should be set to the largest value of frame# observed among all tracks
      if ( tracks_input_RAW(ti).seqOfEvents(end,1) > Nframes )
        Nframes = tracks_input_RAW(ti).seqOfEvents(end,1);
      end
      
      Nevents_RAW(ti) = size( tracks_input_RAW(ti).seqOfEvents, 1 );
      T_RAW(ti)       = tracks_input_RAW(ti).seqOfEvents(Nevents_RAW(ti),1) -  tracks_input_RAW(ti).seqOfEvents(1,1) + 1;
  
      if( size(tracks_input_RAW(ti).tracksCoordAmpCG, 2)/8 ~= (tracks_input_RAW(ti).seqOfEvents(end,1)- tracks_input_RAW(ti).seqOfEvents(1,1) +1 ) )
         disp("inconsistent sizing in Tracks --unclear")
         return
      end
      
      if ( T_RAW(ti) <= 1 )      
          is_ephemeral(ti) = true;
          % if the track only contains a single frame, then we are not interested in it.
      else
          is_ephemeral(ti) = false;
      end
  end
  
  % FROM HERE ON WE ONLY WORK WITH THE non-ephemeral input data  
  
  tracksdat    =  tracks_input_RAW( ~is_ephemeral );
 
  % (nothing with _RAW should be touched after this point )
  
  end