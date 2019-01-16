function [ diffconst_vals_Plist, dpos_Plist, dpos_all, D_obs_all ] =  get_diffdat ( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes ,R )
% diffconst_vals is a data structure (list) containing average diffusion parameters for each state
% dpos_Plist carries the data cataloguing the actual change in position of
% particles.
% D_obs_all (for the histogram) is an array of diffusion constants
% calculated for each track, regardless of state or duration. This is done
% twice (using the old, flawed, method from the PNAS 2013 paper, and the new, correct, way.)

% ----- initialize -----------------

% "Plist" implies a list for [1] monomers, [2] dimers, [3] trimers, etc...
diffconst_vals_Plist  = [];
dpos_Plist            = [];

% also track "all" tracks irrespective of polymer state. dndnp1 is the
% self-correlated movement over 3 successive frames
dpos_all.dndnp1       = [];
dpos_all.dx           = [];
dpos_all.dy           = [];

D_obs_all.oldmethod_slope        = [];
D_obs_all.newmethod_succdx       = [];


%% First collect diffconst_vals for each polymer state

Num_comp_tracks = length( trackdat_xyl);

for S= 1:max_state

    [ dpos_Plist{S}.dx, dpos_Plist{S}.dy ] = get_dxy_list(  state_matrices_allti, trackdat_xyl, Nframes, S);
    MSD_ave               = mean( ( dpos_Plist{S}.dx.^2 + dpos_Plist{S}.dy.^2 ) );

    dxndxnp1_list         = get_dxnnp1_list( state_matrices_allti, trackdat_xyl, Nframes, S);
    dxndxnp1_ave          = mean( dxndxnp1_list );

    % this is the diffusion constant, using the formula from Eq. 14 (page 022726-7
    % from Vestergaard et al, Phys. Rev. E. 89, (2014)
    diffconst_vals_Plist{S}.D       = (MSD_ave/(4*dt)) + dxndxnp1_ave/(dt);
    diffconst_vals_Plist{S}.sigma2  = R*(MSD_ave) + (2*R-1)*dxndxnp1_ave;

    SNR = sqrt( (diffconst_vals_Plist{S}.D + (dxndxnp1_ave/dt))* dt / diffconst_vals_Plist{S}.sigma2 );

    if( SNR < 1  )
        disp( strcat("WARNING: SNR =", num2str( SNR ), " for state ",num2str( S )))
        % disp( "This implies failure to meet criteria of free-diffusion on a non-fluctuating membrane" )
    end

end

% ---------------------------------------------
%-------- now get the overall (stateless) diff constant distribution
%-------- (i.e. calculated for each subtrack)

% minimum frame# length for calculation of D using PNAS method.
Tmin_oldDcalc = 0.1 * Nframes;

for ti = 1:Num_comp_tracks

   % ----------- First collect dx^2 data -----------------------

   mask          =  create_mask( state_matrices_allti{ti}, Nframes, 0, 1 );
   % ----------- the 0 means take any non-NaN state                 ^

   if( ~all( size(mask) == [size( trackdat_xyl(ti).dx, 1), (Nframes-1) ] ) )
      disp("ERROR: mismatched frame length in calc_d2")
      return
   end

   dx2_mat       = mask.*( trackdat_xyl(ti).dx .* trackdat_xyl(ti).dx );
   dy2_mat       = mask.*( trackdat_xyl(ti).dy .* trackdat_xyl(ti).dy );

   d.x            = trackdat_xyl(ti).dx( mask ==1);
   d.y            = trackdat_xyl(ti).dy( mask ==1);
   if( size( d.x,1) ~=1  )
       d.x = transpose( d.x );
       d.y = transpose( d.y );
   end

   dpos_all.dx  = [dpos_all.dx, d.x];
   dpos_all.dy  = [dpos_all.dy, d.y];
   % ----- these will constitute the x-y scatter plot ------------

   % ---- Now collect correlated two-frame sequences -----
   % ---- i.e. dx_n * dx_{n+1}
   mask2             =  create_mask( state_matrices_allti{ti}, Nframes, 0, 2 );

   if( size(mask2,2) ~= Nframes - 2 || size(mask2,1) ~= size( trackdat_xyl(ti).dx, 1) )
      disp("ERROR: mismatched frame length in calc_dnnp1" )
      return
   end

   dxnnp1_mat       = mask2 .*( trackdat_xyl(ti).dx(:,1:Nframes-2) .* trackdat_xyl(ti).dx(:,2:Nframes-1) );
   dynnp1_mat       = mask2 .*( trackdat_xyl(ti).dy(:,1:Nframes-2) .* trackdat_xyl(ti).dy(:,2:Nframes-1) );

   temp   = dxnnp1_mat(mask2 == 1) + dynnp1_mat(mask2 == 1);
   % if there's only a single subtrack, then the array gets returned with a different orientation.
   % Make sure it's a horizonal array:
   if( size( temp,1) ~=1  )
       temp = transpose( temp );
   end

   dpos_all.dndnp1 = [ dpos_all.dndnp1, temp ];

   % ----------- Collect observations of D for all polymer types (D_obs_all)  -----------
   % (will be used for histogram)
   % --- First the new method (based on successive dx values:

   clear MSD;
   Nsubtracks(ti) =  size( dx2_mat, 1);
   for  sti = 1:Nsubtracks(ti)
      temp     = mean(dx2_mat(sti, mask(sti,:)==1 ) + dy2_mat(sti, mask(sti,:)==1 ) );
      if( size(temp,2))
         MSD(sti) = temp;
      else
         MSD(sti) = NaN;
      end
   end

   D_obs_all.newmethod_succdx = [ D_obs_all.newmethod_succdx , (MSD(~isnan(MSD))/(4*dt)) ];
   % this will be our array of "Diffusion constants" for each track observed.
   % irrespective of track length or state.

   clear mask;
   clear MSD;

   % --- Now collect the old (biased) method for comparison:

   for  sti = 1:Nsubtracks(ti)

        valid_frames =  find( ~isnan(trackdat_xyl(ti).Lamp(sti, :) ));

        birth_track_sti = min( valid_frames );
        death_track_sti = max( valid_frames );

        if( death_track_sti - birth_track_sti >= Tmin_oldDcalc && length(valid_frames) > 2 )

            % collect time points after start of track for which
            % luminescance was non-NaN
            valid_frames =  find( ~isnan(trackdat_xyl(ti).Lamp(sti, :) ));
            tvals        = dt* ( valid_frames - birth_track_sti );
            xvals        = trackdat_xyl(ti).xpos(sti, valid_frames ) - trackdat_xyl(ti).xpos(sti,  birth_track_sti );
            yvals        = trackdat_xyl(ti).ypos(sti, valid_frames ) - trackdat_xyl(ti).ypos(sti,  birth_track_sti );

            MSD = ( xvals.^2 + xvals.^2 );
            if ( MSD(1)~= 0 || min(MSD)<0 )
               disp("ERROR: exception in linear regression section.");
            end


            % fit MSD to t values with polynomial of degree n=1
            % C = polyfit(tvals(2:end), MSD(2:end), 1);
            % the slope obtained from this fit is then = 4D
            % x = A\B solves the system of linear equations A*x = B.

            % original paper did not mention any y-intercept.
            slope = tvals(:)\MSD(:);
            % plot (tvals,  xvals.^2 + xvals.^2 , tvals, C(1)*tvals + C(2), tvals, tvals*slope  );

            D_obs_all.oldmethod_slope = [ D_obs_all.oldmethod_slope, slope/4 ];

        end
    end
end

end
