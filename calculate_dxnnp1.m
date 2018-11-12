function [ dnnp1 ] =  calculate_dxnnp1 (statemat, trackmat_xyl, S, Nframes)
% calculate the average squared change in position between adjacent frames 

Num_comp_tracks     = length( trackmat_xyl );

dx_list = []
dy_list = []

% ---------------------------------

for ti = 1:Num_comp_tracks
   % get mask for frames where the track was in the right state in that 
   % frame AND the subsequent one, AND the one after that:
   mask          =  create_mask( statemat{ti}, S, 2 ) 
  
   if( dim(mask,2) ~= Nframes - 2 || dim(mask,1) ~= dim( trackmat_xyl(ti).dx, 1) )
      disp("ERROR: mismatched frame length in calc_dnnp1" )
      return
   end 

   dxnnp1_mat       = mask .*( trackmat_xyl(ti).dx(:,1:Nframes-2) .* trackmat_xyl(ti).dx(:,2:Nframes-1) );
   dynnp1_mat       = mask .*( trackmat_xyl(ti).dy(:,1:Nframes-2) .* trackmat_xyl(ti).dy(:,2:Nframes-1) );

   dxnnp1_list = [ dxnnp1_list , dxnnp1_mat(mask == 1) ];
   dynnp1_list = [ dynnp1_list , dynnp1_mat(mask == 1) ];
   
end

dnnp1 = mean( ( dxnnp1_list + dynnp1_list )/2 );
% Done this way so that if there is a mismatched dimension between
% x and y, the error will be flagged.

