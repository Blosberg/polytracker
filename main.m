% polytracker Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
%
% A software package for processing data from TIRF fluorescance
% videos. Running this program requires that the relevant parameters be set
% as needed below, and that the data structure "TracksFinal" is already
% loaded into the workspace (see documentation at the bottom of this file).


%% ==========================================================================

% ---- DEFINE VARIABLES RELEVANT TO THIS DATASET

dt         = 0.02;       %--- time spacing between frames
px_spacing = 0.106941;   %--- pixel spacing (assuming tracksFinal stores position
                         %    coordinates in units of pixels, this factor converts
                         %    the spatial dimension into micrometers.
R          = 1/6;        %--- the motion blur constant, as defined in Vestergaard et
                         %    al, Phys. Rev. E. 89, (2014). Assuming the camera
                         %    shutter is left on continuously

dim1      = 246;         % length x
dim2      = 183;         % width  (of the field of view, in pixels).

Area       = dim1*dim2*(px_spacing^2);     % This is the cross-sectional area (in um^2) of the field of view
                                           % the number of tracks will be divided by this
                                           % quantity to give you the density per unit area.

Label      = "ds524"; %--- Some descriptive name for your dataset.

Nbin       = 300;     %--- Resolution (number of bins) for your histograms.

% ==================================================================


% ---- GET THE COMPOUND TRACK LENGTHS:
comptrack_lengths = get_comptrack_lengths( tracksFinal );


% Now run the major part of the script:
polydat = polytrack( tracksFinal, Label, dt, px_spacing, R, Area, Nbin);

% ===================================================
% The following documentation kept for reference; taken from the comments of plotComptrack:
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
