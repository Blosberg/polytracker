% polytracker ---a software package for processing data from TIRF fluorescance 
% videos
% Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>

% ---- master file to run the polytracker software. You should be able 
% to hit ctrl+Enter through the commands here after loading your tracksFinal data
% structure into the workspace 

%% ==========================================================================
%% 
% ---- DEFINE VARIABLES RELEVANT TO THIS DATASET

dt         = 0.04;       %--- time spacing between frames
px_spacing = 0.01066;    %--- pixel spacing (assuming tracksFinal stores position
                         %    coordinates in units of pixels, this factor converts 
                         %    the spatial dimension into micrometers.
R          = 1/6;        %--- the motion blur constant, as defined in Vestergaard et 
                         %    al, Phys. Rev. E. 89, (2014). Assuming the camera  
                         %    shutter is left on continuously      

Label      = "DATANAME"; %--- Some descriptive name for your dataset.

Nbin       = 100;        %--- Resolution (number of bins) for your histograms. 

%% ==========================================================================
%% Now run the script: 

[ life, dens, lumen, Diffdat_p, D_observations] = polytrack( tracksFinal, Label, dt, px_spacing, R, Nbin);

