function [ diffconst_vals, dout_Plist, dout_all, dndnp1, D_observations ] =  get_diffdat ( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes ,R )
% diffconst_vals is a data structure (list) containing average diffusion parameters for each state
% d2out_Plist is an x-y array of all the squared changes in position for each state --likewise for _all
% except binned together into a single list.
% dndnp1 is an array of correlations between adjacent steps. Checking for drift.
% D_observations is an array of diffusion constants calculated for each track, regardless of state or duration.

% ----- initialize -----------------
D_observations  = [];
dout_Plist      = [];
dout_all.x      = [];
dout_all.y      = [];
dndnp1          = [];

% -------------------------------------------------
% ---- First collect diffconst_vals for each polymer state

Num_comp_tracks = length( trackdat_xyl);

for S= 1:max_state

    [ dout_Plist{S}.dx, dout_Plist{S}.dy ] = get_dxy_list(  state_matrices_allti, trackdat_xyl, Nframes, S);
    MSD_ave               = mean( ( dout_Plist{S}.dx.^2 + dout_Plist{S}.dy.^2 ) );

    dxndxnp1_list         = get_dxnnp1_list( state_matrices_allti, trackdat_xyl, Nframes, S);
    dxndxnp1_ave          = mean( dxndxnp1_list );

    % this is the diffusion constant, using the formula from Eq. 14 (page 022726-7
    % from Vestergaard et al, Phys. Rev. E. 89, (2014)
    diffconst_vals{S}.D       = (MSD_ave/(4*dt));
    diffconst_vals{S}.sigma2  = R*(MSD_ave) + (2*R-1)*dxndxnp1_ave;

    SNR = sqrt( (diffconst_vals{S}.D + (dxndxnp1_ave/dt))* dt / diffconst_vals{S}.sigma2 );

    if( SNR < 1  )
        disp( strcat("WARNING: SNR =", num2str( SNR ), " for state ",num2str( S )))
        disp( "This implies failure to meet criteria of free-diffusion on a non-fluctuating membrane" )
    end

end

% ---------------------------------------------
%-------- now get the overall (stateless) diff constant distribution
%-------- (i.e. calculated for each subtrack)

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

   D_observations = [ D_observations , (MSD(~isnan(MSD))/(4*dt)) ];
   % this will be our array of "Diffusion constants" for each track observed.
   % irrespective of track length or state.

   d.x            = trackdat_xyl(ti).dx( mask ==1);
   d.y            = trackdat_xyl(ti).dy( mask ==1);
   if( size( d.x,1) ~=1  )
       d.x = transpose( d.x );
       d.y = transpose( d.y );
   end


   dout_all.x  = [dout_all.x, d.x];
   dout_all.y  = [dout_all.y, d.y];
   % ----- these will constitute the x-y scatter plot ------------

   clear mask;
   % ---- now consider correlations with two-frame sequences -----
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

   dndnp1 = [ dndnp1, temp ];

end

end
