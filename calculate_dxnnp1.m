function [ dnnp1 ] =  calculate_dxnnp1 (state_matrices_allti, trackdat_xyl, Nframes, S)
% collect change in position between two adjacent frames

Num_comp_tracks     = length( trackdat_xyl );

dxnnp1_list  = [];
dynnp1_list  = [];

% ---------------------------------

for ti = 1:Num_comp_tracks
   % get mask for frames where the track was in the right state in that
   % frame AND the subsequent one, AND the one after that:
   mask          =  create_mask( state_matrices_allti{ti}, Nframes, S, 2 );

   if( size(mask,2) ~= Nframes - 2 || size(mask,1) ~= size( trackdat_xyl(ti).dx, 1) )
      disp("ERROR: mismatched frame length in calc_dnnp1" )
      return
   end

   dxnnp1_mat       = mask .*( trackdat_xyl(ti).dx(:,1:Nframes-2) .* trackdat_xyl(ti).dx(:,2:Nframes-1) );
   dynnp1_mat       = mask .*( trackdat_xyl(ti).dy(:,1:Nframes-2) .* trackdat_xyl(ti).dy(:,2:Nframes-1) );

   dxnnp1_vec = dxnnp1_mat(mask == 1);
   dynnp1_vec = dynnp1_mat(mask == 1);

   % if there's only a single subtrack, then the array gets returned with a different orientation.
   % Make sure it's a horizonal array:
   if( size(dxnnp1_vec,1) ~=1  )
       dxnnp1_vec = transpose( dxnnp1_vec );
       dynnp1_vec = transpose( dynnp1_vec );
   end

  dxnnp1_list = [ dxnnp1_list , dxnnp1_vec ];
  dynnp1_list = [ dynnp1_list , dynnp1_vec ];


end

dnnp1 = mean( ( dxnnp1_list + dynnp1_list )/2 );
% if there is any mismatched dimension between x and y,
% then there is a larger-scale problem that we want flagged.

end
