function [ Diffdat, dout, dndnp1, D_observations ] =  get_diffdat ( state_matrices_allti, trackdat_xyl, max_state, dt, Nframes ,R )
% Diffdat is a data structure (list) containing diffusion parameters for each state
% dout is an x-y array of all the changes in position. Checking for drift
% dndnp1 is an array of correlations between adjacent steps. Checking for drift.
% D_observations is an array of diffusion constants calculated for each track, regardless of state or duration.

% ----- initialize -----------------
D_observations  = [];
dout.x          = [];
dout.y          = [];
dndnp1          = [];

% -------------------------------------------------
% ---- First collect Diffdat for each polymer state

Num_comp_tracks = length( trackdat_xyl);

for S= 1:max_state

    dx2               = calculate_dx2(  state_matrices_allti, trackdat_xyl, Nframes, S);
    dx2_ave           = mean(dx2);

    dxndxnp1          = calculate_dxnnp1( state_matrices_allti, trackdat_xyl, Nframes, S);
    dxndxnp1_ave      = mean(dxndxnp1);

    % this is the diffusion constant, using the formula from Eq. 14 (page 022726-7
    % from Vestergaard et al, Phys. Rev. E. 89, (2014)
    Diffdat{S}.D       = (dx2_ave/2*(dt)) + (dxndxnp1_ave/dt);
    Diffdat{S}.sigma2  = R*(dx2_ave) + (2*R-1)*dxndxnp1_ave;

    if( Diffdat{S}.sigma2 <1  )
        disp( strcat("WARNING: sigma^2 =", num2str(Diffdat{S}.sigma2 ), " for state ",num2str( S ) ) )
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

   
   dout.x  = [dout.x, d.x];
   dout.y  = [dout.y, d.y];
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