function [ dx_list, dy_list ] =  get_dxy_list (state_matrices_allti, trackdat_xyl, Nframes, S)
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
      disp("ERROR: mismatched frame length in get_d2_list")
      return
   end

   dx_mat       = mask.*( trackdat_xyl(ti).dx );
   dy_mat       = mask.*( trackdat_xyl(ti).dy );

   dx_vec = dx_mat(mask == 1);
   dy_vec = dy_mat(mask == 1);

   % if there's only a single subtrack, then the array gets returned with a different orientation.
   % Make sure it's a horizonal array:
   if( size(dx_vec,1) ~=1  )
       dx_vec = transpose( dx_vec );
       dy_vec = transpose( dy_vec );
   end

   dx_list = [ dx_list , dx_vec ];
   dy_list = [ dy_list , dy_vec ];

end

% d2_list = (dx_list+dy_list)/2 ;
% if there is any mismatched dimension between x and y,
% then there is a larger-scale problem that we want flagged.

end
