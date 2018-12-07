% polytracker ---a software package for processing data from TIRF fluorescance 
% videos
% Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>

% last updated Dec. 7 2018
% ---- master file to run the polytracker software. You should be able 
% to hit ctrl+Enter through the commands here after loading your tracksFinal data
% structure into the workspace 

%% ==========================================================================
%% 
% ---- DEFINE VARIABLES RELEVANT TO THIS DATASET

dt         = 0.04;       %--- time spacing between frames
px_spacing = 0.106941;   %--- pixel spacing (assuming tracksFinal stores position
                         %    coordinates in units of pixels, this factor converts 
                         %    the spatial dimension into micrometers.
R          = 1/6;        %--- the motion blur constant, as defined in Vestergaard et 
                         %    al, Phys. Rev. E. 89, (2014). Assuming the camera  
                         %    shutter is left on continuously      

Area       = 144;       % this is the cross-sectional area of the field of view
                         % the number of tracks will be divided by this
                         % quantity to give you the density per unit area.

Label      = "Before-DAMGO"; %--- Some descriptive name for your dataset.

Nbin       = 200;        %--- Resolution (number of bins) for your histograms. 

%% ==========================================================================
%% Now run the script: 

[ lifetime_list, density, lumen_list, dout_Plist, Diffconst_vals, D_observations] = polytrack( tracksFinal, Label, dt, px_spacing, R, Area, Nbin);

% -----------------------------
% and plot the results:

max_state = length(lifetime_list);
all_dimerizations = [];
for s = 2:max_state
   mean_lifetime = mean( lifetime_list{s}) ;
   all_dimerizations = [ all_dimerizations, lifetime_list{s}];
end

figure(1);
hist(all_dimerizations, 50);
xlabel(" lifetime [s]")
ylabel("frequency")
title( strcat('Merger liftime distribution; dataset: ', Label) )
