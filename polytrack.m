function polydat =  polytrack ( tracks_input_RAW, Label, dt, px_spacing, R, Area, Nbin)
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


% Previous code instances were interested in the lifetime of various Polymer states.
% this is left here in case we want this functionality again later
% polydat.lifetime_list = get_state_lifetimes ( tracks_input, state_matrices_allti, max_state, dt, Nframes );

[ polydat.bond_time.on, polydat.bond_time.off ] = get_bond_times (  tracks_input, state_matrices_allti, max_state, dt, Nframes  );


figure(1);
subplot(2,1,1);
hist( polydat.bond_time.on, 50);
xlabel(" On [s]")
ylabel("frequency")
title( strcat('On/Off window liftime distribution; dataset: ', Label) )


subplot(2,1,2);
hist( polydat.bond_time.off, 50);
xlabel(" Off [s]")
ylabel("frequency")

%%  ===================================================
% Collect density
polydat.density = get_particle_density( state_matrices_allti, max_state, Nframes, Area );


figure(2)
tvals = dt*linspace(1,Nframes,Nframes);

plot(tvals, sum(polydat.density,1) )
xlabel("time")
ylabel("density of distinct tracks")
title(  strcat('number of distinct tracks; dataset:', Label) )

figure(3)
plot(tvals, sum( (1:max_state)'.* polydat.density, 1) )
xlabel("time")
ylabel("cumulative density of polymers")
title( strcat('Number of proteins total in frame (summed over all states); dataset: ', Label) )

figure(4)
plot(tvals, polydat.density(1,:) ./ sum(polydat.density,1) )
xlabel("time")
ylabel("Fraction of monomers")
title( strcat('Fraction of monomers/tracks; dataset: ', Label) )


% ===================================================
% --- Tabulate luminescance by state
% same convention as above:


polydat.lumen_Plist = get_lumen_list ( tracks_input, state_matrices_allti, trackdat_xyl, max_state, Nframes );

for s = 1:max_state
   mean_lumen(s) = mean( polydat.lumen_Plist{s} );
   lumen_stdv(s) = std(  polydat.lumen_Plist{s} );
end

figure(5)
plot( mean_lumen - lumen_stdv,"--r"); hold on;
plot( mean_lumen + lumen_stdv,"--r"); hold on;
plot( mean_lumen,'r', 'LineWidth', 4  );
axis([0 max_state min(mean_lumen - lumen_stdv)  max(mean_lumen + lumen_stdv) ])
xlabel("state")
ylabel("mean illumination intensity +/- std. dev")
title( strcat('Mean illumination by state (1=monomer, 2=dimer, etc.); dataset: ', Label) )


figure(6)
subplot(2,1,1);
hist( polydat.lumen_Plist{1}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of monomers; dataset: ', Label) )
% xlim([0, 0.01]);

subplot(2,1,2);
hist( polydat.lumen_Plist{2}, Nbin)
xlabel("Intensity")
ylabel("Freq")
title( strcat('spectral distribution of dimers; dataset: ', Label) )
% xlim([0, 0.01]);



% ===================================================
% get diffusion constants:
% same convention as above:


[ polydat.diffus.Dvals_Plist, polydat.diffus.dpos_Plist, dpos_all, polydat.diffus.D_observations ]  = get_diffdat( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes, R  );

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
hist( polydat.diffus.D_observations.oldmethod_slope, 2*Nbin);
xlabel("D");
ylabel("Freq");
title( strcat('Diffusion constant calculation as in PNAS 2013: ', Label) );
% xlim(Dval_lims);

subplot(2,1,2);
hist( polydat.diffus.D_observations.newmethod_succdx, 2*Nbin);
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
