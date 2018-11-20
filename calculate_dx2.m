function [ d2 ] =  calculate_dx2 (state_matrices_allti, trackdat_xyl, Nframes, S)
% collect squared change in position between adjacent frames

Num_comp_tracks     = length( trackdat_xyl );

% initialize results
dx_list = [];
dy_list = [];
% ---------------------------------

for ti = 1:Num_comp_tracks

   % get mask for frames where the track was in the right state in that
   % frame AND the subsequent one:
   mask          =  create_mask( state_matrices_allti{ti}, Nframes, S, 1 );

   if( ~all( size(mask) == [size( trackdat_xyl(ti).dx, 1), (Nframes-1) ] ) )
      disp("ERROR: mismatched frame length in calc_d2")
      return
   end

   dx2_mat       = mask.*( trackdat_xyl(ti).dx .* trackdat_xyl(ti).dx );
   dy2_mat       = mask.*( trackdat_xyl(ti).dy .* trackdat_xyl(ti).dy );

   dx2_vec = dx2_mat(mask == 1);
   dy2_vec = dy2_mat(mask == 1);

   % if there's only a single subtrack, then the array gets returned with a different orientation.
   % Make sure it's a horizonal array:
   if( size(dx2_vec,1) ~=1  )
       dx2_vec = transpose( dx2_vec );
       dy2_vec = transpose( dy2_vec );
   end

   dx_list = [ dx_list , dx2_vec ];
   dy_list = [ dy_list , dy2_vec ];

end

d2 = mean( (dx_list+dy_list)/2  );
% if there is any mismatched dimension between x and y,
% then there is a larger-scale problem that we want flagged.

end

