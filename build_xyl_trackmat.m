function [ trackmat_xyl ] =  build_xyl_trackmat ( tracks_input, Nframes )

Num_comp_tracks = length( tracks_input);

for ti = 1:Num_comp_tracks

    Nsubtracks(ti)  =  size( tracks_input(ti).tracksCoordAmpCG, 1);

    % get the offset from which this composite track starts and finishes.
    comptrack_birth     = min( tracks_input(ti).seqOfEvents(:,1) );

    comptrack_death_lag = Nframes - max( tracks_input(ti).seqOfEvents(:,1) );
    % the number of frames that should be appended to the end of this comptrack to make it to the end of the movie.

    % Loop through each independent track (of length T frames) and figure
    % out how many sub-tracks need to be considered for each case:
    for sti = 1:Nsubtracks( ti )
       % Within each sub-track, extract the position and change values
       % segment in which the particle could not be located on either side
       % of the window

       % get birth and death time points for each subtrack _within_ the composite window
       obituary = tracks_input(ti).seqOfEvents( tracks_input(ti).seqOfEvents(:,3)==sti, 1 );

       xdat_compwindow = tracks_input(ti).tracksCoordAmpCG(sti,1:8:end);
       ydat_compwindow = tracks_input(ti).tracksCoordAmpCG(sti,2:8:end);
       Adat_compwindow = tracks_input(ti).tracksCoordAmpCG(sti,4:8:end);

       trackmat_xyl(ti).xpos(sti,:) = [ repmat(NaN,[1 , (comptrack_birth-1)]), xdat_compwindow, repmat(NaN,[1 , comptrack_death_lag]) ];
       trackmat_xyl(ti).ypos(sti,:) = [ repmat(NaN,[1 , (comptrack_birth-1)]), ydat_compwindow, repmat(NaN,[1 , comptrack_death_lag]) ];
       trackmat_xyl(ti).Lamp(sti,:) = [ repmat(NaN,[1 , (comptrack_birth-1)]), Adat_compwindow, repmat(NaN,[1 , comptrack_death_lag]) ];

    end % terminate forloop through subtracks

    if(  ~all(  size( trackmat_xyl(ti).xpos ) == [ Nsubtracks(ti), Nframes] ) ||  ~all(  size( trackmat_xyl(ti).ypos ) == [ Nsubtracks(ti), Nframes] ) )
       disp("unexpected dimension in trackmat_xyl.")
       return
    end

    trackmat_xyl(ti).dx  = diff ( trackmat_xyl(ti).xpos, 1, 2 );
    trackmat_xyl(ti).dy  = diff ( trackmat_xyl(ti).ypos, 1, 2 );

    if(  ~all (  size( trackmat_xyl(ti).Lamp ) ==  [ Nsubtracks(ti), Nframes] ) )
        disp("unexpected dimension in trackmat_xyl.Lamp.")
        return
    end

end % terminate forloop through comptracks


end
