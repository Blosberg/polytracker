# polytracker
This matlab package builds upon fluorescence-tracking data from other packages such as u-track (see https://github.com/mcianfrocco/Matlab/), and streamlines the processing of various relevent observables.

Users should open the file "main.m" and edit the variables defined there, as needed, to correspond to their input data. E.g. frame spacing (dt), field-of-view ares, etc., and some "Label" to assign to your data in the resulting plots.

A data structure like "TracksFinal" should also be loaded into the workspace. From that point, the main script can simply be run (e.g. hit "ctrl+Enter" from the matlab IDE, and the package will run; various observables are plotted automatically.


The main novelty of the package is inferring the state of tracks based on their merger/splitting events. As such, most observables are output with the suffix "Plist" --this implies observations specific to [1] monomers, [2] dimers, [3] trimers, etc. 

N.B. with higher-order oligomerizations, the ability of polytracker to accurately determine the state via events diminishes. A useful indicator for quality-control purposes is plotted automatically comparing the assigned state to average luminescance intensity.  The slope of this line will be biased downward due to preferential photo-bleaching of higher-states (e.g. dimer lumuniscence will be less than double that of monomers), but the general upward trend should be robust. Rather than merely comparing dimers and monomers on average, however, polytracker assigns state values to individual tracks. 

Diffusion calculations are performed using the theory outlined in Vestergaard et al, Phys. Rev. E. 89, (2014), (n.b. Eq. 14 (page 022726-7)) to avoid the known biases of long-term trendlines in MSD vs. time. scatter plots and histograms are produced to illustrate trends in particle movement, and to check for persistent (non-Brownian) movement
Histograms of diffusion data across all polymer states are produced, and individual diffusion values per polymer state (i.e. as a _Plist construct) are under development.


