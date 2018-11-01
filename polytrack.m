% last edited from Boontower 17:06, 25.10.18
% Eventually turn this into a function
% function [velocity_batch] = analyze_track(tracks_input_RAW, dt )

%===== remove this when we make this a funciton =====
% dat_in = "E:\NikonTIRF\04-10-18\beta1\141\TrackingPackage\tracks\Channel_1_tracking_result"
% tracksoftware = "C:\u-track\software\plotCompTrack.m"
% load(dat_in)
% clear *

load("workspace_just_tracksFinal.m")

tracks_input_RAW = tracksFinal; 
dt     = 0.04
% ===================================================
% --- FILTER OUT EPHEMERAL TRACKS:
% select only the tracks that have length greater than 1 frame.
% ne == "non-ephemeral" --i.e. of finite duration.
% can probably simplify this into a one-liner
% something like:
% is_ephemeral = length( tracks_input.tracksCoordAmpCG ) > 8
% @@@ try to remove the T

Num_indep_tracks_RAW = length(tracks_input_RAW);
Nframes = 1

for ti = 1:Num_indep_tracks_RAW
    
    # Nframes should be set to the largest value of frame# observed among all tracks
    if ( tracksFinal(ti).seqOfEvents(end,1) > Nframes )
      Nframes = tracksFinal(ti).seqOfEvents(end,1);
    end
    
    Nevents_RAW(ti) = size( tracks_input_RAW(ti).seqOfEvents, 1 );
    T_RAW(ti)       = tracks_input_RAW(ti).seqOfEvents(Nevents_RAW(ti),1) -  tracks_input_RAW(ti).seqOfEvents(1,1) + 1;

    if( size(tracks_input_RAW(ti).tracksCoordAmpCG, 2)/8 ~= (tracksFinal(ti).seqOfEvents(end,1)- tracksFinal(ti).seqOfEvents(1,1) +1 ) )
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

tracks_input     =  tracks_input_RAW( ~is_ephemeral );
T                =  T_RAW( ~is_ephemeral );      % -- Number of Frames for that compound track
Nevents          =  Nevents_RAW( ~is_ephemeral ); % -- Number of events in that compound track

Num_indep_tracks = length( tracks_input);

% (nothing with _RAW should be touched after this point )

% ===================================================
% --- get xy- and dxdy- data for each subtrack

for ti = 1:Num_indep_tracks

    Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);   

    % Loop through each independent track (of length T frames) and figure
    % out how many sub-tracks need to be considered for each case:
    for sti = 1:Nsubtracks( ti )
       % Within each sub-track, extract the position and change values
       % the dxdy array will be of size T(ti)-1 and have values NA for any
       % segment in which the particle could not be located on either side
       % of the window   
    
       trackdat_xy(ti).subtrack(sti).xypos(1,:) = tracks_input(ti).tracksCoordAmpCG(sti,1:8:end);
       trackdat_xy(ti).subtrack(sti).xypos(2,:) = tracks_input(ti).tracksCoordAmpCG(sti,2:8:end);
      
       trackdat_xy(ti).subtrack(sti).dxdy(1,:)  = diff ( trackdat_xy(ti).subtrack(sti).xypos(1,:) );
       trackdat_xy(ti).subtrack(sti).dxdy(2,:)  = diff ( trackdat_xy(ti).subtrack(sti).xypos(2,:) );
    end

end

% ===================================================
% --- Assign polymer state for each time point along each track:
% NaN = non-existence, 1 = monomer, 2 = dimer, 3 = trimer,  etc...

max_state = 1
for ti = 1:Num_indep_tracks % go through all non-ephemeral independent track sets.
   
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
                state{ti}( tracks_input(ti).seqOfEvents(evi,3) , (tracks_input(ti).seqOfEvents(evi,1)): Nframes ) = NaN ;
                %                                           ^ col. 3 =subtrack that just died.
            end
           
            % ---- death through merger ?
            if ( ~isnan( tracks_input(ti).seqOfEvents(evi,4) )  )
                % if so, then increment state of the other track affected by this
                % point on to the end of the full track.
                temp_inc = state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), tracks_input(ti).seqOfEvents(evi,1) : Nframes );
                           state{ti}( tracks_input(ti).seqOfEvents(evi,4 ), tracks_input(ti).seqOfEvents(evi,1) : Nframes ) = temp_inc + 1 ;
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
    end
   
    temp = max( max( state{ti} ) );
    if (temp > max_state)
       max_state = temp;
    end

end % --- finished for-loop over ti through all non-ephemeral independent track sets.

%% 
% ===================================================
% get list of arrays of lifetimes observed for each polymer state  --> MS
lifetime_list = get_state_lifetimes ( state, tracks_input, max_state, Nframes );

% ===================================================

for ti = 1:Num_indep_tracks
     
      mask       = create_mask( state_matrix{ti}, S, 0 )
      maskA      = resize( mask, 1, ... ) @@@

% simple way to convert matrix to array:
% yourvector = yourmatrix(:); @@@ test this.

      deriv_mask = get_state_mask( state_matrix{ti}, S, 1 )
      deriv_maskA= resize( deriv_maskA, S, 0 )

      Lumen_temp = resize( mask.*Lumen{ti} ) @@@ grab lumen data above
      Polydat{S}.Lumen = [ Polydat{S}.Lumen, Lumen_temp[ maskA != 0 ] ]
    
      dx_temp       = resize( mask.*dx{ti} ) @@@ grab dx data above
      Polydat{S}.dx = [ Polydat{S}.Lumen, Lumen_temp[ maskA != 0 ] ]

      dy_temp       = resize( mask.*dy{ti} ) @@@ grab dy data above
      Polydat{S}.dy = [ Polydat{S}.Lumen, Lumen_temp[ maskA != 0 ] ]

end


% ===================================================
% get diffusion constants:

for S = 1:MS
  
   Polydat{S}.diff_const = get_diff_const( Polydat{S}.dx, Polydat{S}.dy,  )

end

% for the diffusion constant, consider these functions:
%  http://tinevez.github.io/msdanalyzer/
% 
%  https://de.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker
% 


% ===================================================
% Take histograms of all the above variable:

% figure(1);
% plot Lumen data for mono, di, tri-mer

% figure(2);
% plot dx,dy data for diffusion const

% figure(3);
% plot histograms of lifetimes for each type of polymer
% histogram(X,nbins)

% ===================================================
% --- build masks



% ===================================================
% --- Tabulate luminescance by state
% Lumen_dat{1} is the list of light intensity observations for polymers in the 1 state
% Lumen_dat{2} "" "" in the 2 state, etc.
% e.g. nameArray = [nameArray, 'Name you want to append'];
% or c = cat(1,a,b)


% ===================================================
% ---

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


