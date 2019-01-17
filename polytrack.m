function [ lifetime_list, Diffconst_vals_Plist, dpos_Plist, lumen_Plist, density, D_observations] =  polytrack ( tracks_input_RAW, Label, dt, px_spacing, R, Area, Nbin)
% take TracksFinal type data structure and dt spacing and some label, and
% perform relevant analysis with a few output plots

%%  ===================================================
% ---  Filter out ephemeral tracks:

[ tracks_input, Nframes ]  = purge_ephemeral_and_disordered( tracks_input_RAW );
% no reference to tracks_input_RAW should ever be made after this point. 
% We now use only tracks_input, where problematic tracks have been removed.

%%  ===================================================
% --- get xy-, dxdy-, and lumen data for each subtrack

trackdat_xyl =  build_xyl_trackdat ( tracks_input, px_spacing, Nframes );

% output arrays in trackdat_xyl include:  .xpos[ition], ypos[ition], 
% and the change thereof( .dx, .dy)
% and luminescance amplitude (.Lamp)

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

density.monomers = (1/Area)* density.monomers;
density.tracks   = (1/Area)* density.tracks;
density.weighted_polymers = (1/Area)* density.weighted_polymers;

figure(2)
tvals = dt*linspace(1,Nframes,Nframes);

plot(tvals, density.tracks)
xlabel("time")
ylabel("density of distinct tracks")
title(  strcat('number of distinct tracks; dataset:', Label) )

figure(3)
plot(tvals, density.weighted_polymers)
xlabel("time")
ylabel("cumulative density of polymers")
title( strcat('Number of proteins total in frame (summed over all states); dataset: ', Label) )

figure(4)
plot(tvals, density.monomers./density.tracks)
xlabel("time")
ylabel("Fraction of monomers")
title( strcat('Fraction of monomers/tracks; dataset: ', Label) )

% ===================================================
% --- Tabulate luminescance by state
% same convention as above:

lumen_Plist = get_lumen_list ( tracks_input, state_matrices_allti, trackdat_xyl, max_state, Nframes );

for s = 1:max_state
   mean_lumen(s) = mean( lumen_Plist{s} );
end

figure(5)
plot( mean_lumen )
xlabel("state")
ylabel("mean illumination intensity")
title( strcat('Mean illumination by state (1=monomer, 2=dimer, etc.); dataset: ', Label) )


figure(6)
subplot(2,1,1);
hist(lumen_Plist{1}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of monomers; dataset: ', Label) )
% xlim([0, 0.01]);

subplot(2,1,2);
hist(lumen_Plist{2}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of dimers; dataset: ', Label) )
% xlim([0, 0.01]);


% ===================================================
% get diffusion constants:
% same convention as above:


[ Diffconst_vals_Plist, dpos_Plist, dpos_all, D_observations ]  = get_diffdat( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes, R  );

% ---------
figure(7);
subplot(2,1,1);
plot(dpos_all.dx, dpos_all.dy, '.');
xlabel("dx");
ylabel("dy");
title( strcat('position change -scatter') );

subplot(2,1,2);
hist( dpos_all.dndnp1, 2*Nbin);
%xlim([-0.05, 0.05]);
xlabel("dx_n * dx_{n+1}");
ylabel("Freq");
title( strcat('Two-frame drift correlation: ', Label) );

% ---------
Dval_lims = [0, 0.2];

figure(8);
subplot(2,1,1);
hist( D_observations.oldmethod_slope, 2*Nbin);
xlabel("D");
ylabel("Freq");
title( strcat('Diffusion constant calculation as in PNAS 2013: ', Label) );
% xlim(Dval_lims);

subplot(2,1,2);
hist( D_observations.newmethod_succdx, 2*Nbin);
xlabel("D");
ylabel("Freq");
title( strcat('Diffusion constant calculation new method: ', Label) );
% xlim(Dval_lims);

% for the diffusion constant, consider these functions:
% http://tinevez.github.io/msdanalyzer/
%
% https://de.mathworks.com/matlabcentral/fileexchange/34040-simple-tracker
%

end
