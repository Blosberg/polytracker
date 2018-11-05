% Eventually turn this into a function
% function [ trackinfo ] =  polytrack ( tracks_input, Label )

% ===================================================
% dat_in = "E:\NikonTIRF\04-10-18\beta1\141\TrackingPackage\tracks\Channel_1_tracking_result"
% tracksoftware = "C:\u-track\software\plotCompTrack.m"
% load(dat_in)
% clear *


load("workspace_just_tracksFinal.m")
tracks_input_RAW = tracksFinal; 
dt     = 0.04
% ===================================================


% --- FILTER OUT EPHEMERAL TRACKS:
[ tracks_input, Nframes ]  = purge_ephemeral( tracks_input_RAW );

% ===================================================
% --- get xy-, dxdy-, and lumen data for each subtrack

trackdat_xyl =  build_xyl_trackmat ( tracks_input, Nframes );

% ===================================================
% --- Assign polymer state for each time point along each track:
% NaN = non-existence, 1 = monomer, 2 = dimer, 3 = trimer,  etc...

[ state, max_state ] =  build_state_matrix ( tracks_input, trackdat_xyl,  Nframes );
 
% ===================================================
% get list of arrays of lifetimes observed for each polymer state  --> MS
% Lifetime_list{1} is the list of lifetime observations for polymers in the 1 state
% Lifetime_list{2} "" "" in the 2 state, etc.

lifetime_list = get_state_lifetimes ( tracks_input, state, max_state, Nframes );
all_lives_matter = []

for s = 2:max_state
   mean_lifetime = mean( lifetime_list{s}) ;
   all_lives_matter = [ all_lives_matter, lifetime_list{s}];
end

figure(1)
hist(all_lives_matter, 50)
xlabel("# frames")
ylabel("frequency")
title('Merger liftime distribution')

% ===================================================
% Collect density 
[ tracknum_density, weighted_polymer_dens ] = get_particle_density( state, Nframes )
figure(2)
plot(tracknum_density)
xlabel("Frame index")
ylabel("# of distinct tracks")
title('number of live, distinct tracks vs. time. \n(ignoring state)')

figure(3)
plot(weighted_polymer_dens)
xlabel("Frame index")
ylabel("cumulative polymer presence")
title('Number of monomers in frame (summed over all states)')


% ===================================================
% --- Tabulate luminescance by state
% same convention as above:
 
lumen_list = get_lumen_list ( tracks_input, state, trackdat_xyl, max_state, Nframes );

for s = 1:max_state
   mean_lumen = mean( lumen_list{s})  
end


% ===================================================
% get diffusion constants:
% same convention as above:

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

% figure(2);
% plot dx,dy data for diffusion const

% figure(3);
% plot histograms of lifetimes for each type of polymer
% histogram(X,nbins)

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


