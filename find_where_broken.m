function [ broken_track, breaking_event ] =  find_where_broken ( SoE_statemat_in )
% Input: statemat. Find a place where states become <= 0
% should be the track where it happens _last_ 

% [ frame#,  1/2=birth/death, THIS track # , othertrack split/merge with ]

[ Nevents , Nsubtracks ] = size( SoE_statemat_in );

breaking_event = 0;

% scan through each track looking for untenables
for Tri = 1:Nsubtracks

    % get the indices of non-positive numeric state values for this track
    Trseq = SoE_statemat_in(:, Tri);
    unten = find ( Trseq < 1 ); % TODO @@@ check that NaN < 1 returns F
   
    % if any occure: 
    if ( length(unten) > 1 )
       % then take the first such index
       unten_onset = min( unten );
       % and if it occurs later than the current last one, 
       % set it as our breaking event.
       if ( unten_onset > breaking_event )
          breaking_event = unten_onset;
          broken_track   = Tri;
       end 
    end

end

if ( breaking_event <= 0 )
   disp("ERROR in find_where_broken. Did not find anything broken.")
end

end
