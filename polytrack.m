% take TracksFinal type data structure and dt spacing and some label, and
% output relevant plots

function [ lifetime_list, density, lumen_list, dout_Plist, Diffconst_vals, D_observations] =  polytrack ( tracks_input_RAW, Label, dt, px_spacing, R, Nbin)

% ===================================================
% dat_in = "E:\NikonTIRF\04-10-18\beta1\141\TrackingPackage\tracks\Channel_1_tracking_result"
% tracksoftware = "C:\u-track\software\plotCompTrack.m"
% load(dat_in)
% clear *
% load("workspace_just_tracksFinal.mat")
% tracks_input_RAW = tracksFinal;
% dt     = 0.04;

%%  ===================================================


% --- FILTER OUT EPHEMERAL TRACKS:
[ tracks_input, Nframes ]  = purge_ephemeral( tracks_input_RAW );


%%  ===================================================
% --- get xy-, dxdy-, and lumen data for each subtrack

trackdat_xyl =  build_xyl_trackdat ( tracks_input, px_spacing, Nframes );

%%  ===================================================
% --- Assign polymer state for each time point along each track:
% NaN = non-existence, 1 = monomer, 2 = dimer, 3 = trimer,  etc...

[ state_matrices_allti, max_state ] =  build_state_matrix ( tracks_input, trackdat_xyl,  Nframes );

%%  ===================================================
% get list of arrays of lifetimes observed for each polymer state  --> MS
% Lifetime_list{1} is the list of lifetime observations for polymers in the 1 state
% Lifetime_list{2} "" "" in the 2 state, etc.

lifetime_list = get_state_lifetimes ( tracks_input, state_matrices_allti, max_state, dt, Nframes );


%%  ===================================================
% Collect density
[ density.monomers, density.tracks, density.weighted_polymers ] = get_particle_density( state_matrices_allti, Nframes );
figure(2)
tvals = dt*linspace(1,Nframes,Nframes);

plot(tvals, density.tracks)
xlabel("time")
ylabel("# of distinct tracks")
title(  strcat('number of distinct tracks; dataset:', Label) )

figure(3)
plot(tvals, density.weighted_polymers)
xlabel("time")
ylabel("cumulative polymer presence")
title( strcat('Number of proteins total in frame (summed over all states); dataset: ', Label) )

figure(4)
plot(tvals, density.monomers./density.tracks)
xlabel("time")
ylabel("Fraction of monomers")
title( strcat('Fraction of monomers/tracks; dataset: ', Label) )

% ===================================================
% --- Tabulate luminescance by state
% same convention as above:

lumen_list = get_lumen_list ( tracks_input, state_matrices_allti, trackdat_xyl, max_state, Nframes );

for s = 1:max_state
   mean_lumen(s) = mean( lumen_list{s});
end

figure(5)
plot( mean_lumen )
xlabel("state")
ylabel("mean illumination intensity")
title( strcat('Mean illumination by state (1=monomer, 2=dimer, etc.); dataset: ', Label) )


figure(6)
subplot(2,1,1);
hist(lumen_list{1}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of monomers; dataset: ', Label) )

subplot(2,1,2);
hist(lumen_list{2}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of dimers; dataset: ', Label) )

% ===================================================
% get diffusion constants:
% same convention as above:


[ Diffconst_vals, dout_Plist, dout_all, dndnp1, D_observations ]  = get_diffdat( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes, R  );

figure(7);
subplot(2,1,1);
plot(dout_all.x, dout_all.y, '.');
xlabel("dx");
ylabel("dy");
title( strcat('position change -scatter') );

subplot(2,1,2);
hist(dndnp1, 2*Nbin);
xlim([-0.0005, 0.0005]);
xlabel("dx_n * dx_{n+1}");
ylabel("Freq");
title( strcat('Two-frame drift correlation: ', Label) );

figure(8);
hist(D_observations, 2*Nbin);
xlabel("Observed diffusion constant");
ylabel("Freq");
title( strcat('Diffusion constant calculation (like in PNAS 2013): ', Label) );
xlim([0, 0.0050]);

% for the diffusion constant, consider these functions:
%  http://tinevez.github.io/msdanalyzer/
%
%  https://de.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker
%

end
