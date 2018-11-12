function [ d2 ] =  calculate_dx2 (statemat_list, trackmat_xyl, S, Nframes)
% calculate the average squared change in position between adjacent frames 

Num_comp_tracks     = length( trackmat_xyl );

dx_list = []
dy_list = []

% ---------------------------------

for ti = 1:Num_comp_tracks
   % get mask for frames where the track was in the right state in that 
   % frame AND the subsequent one:
   mask          =  create_mask( statemat_list{ti}, S, 1 );
  
   if( dim(mask,2) ~= Nframes - 1 || dim(mask,1) ~= dim( trackmat_xyl(ti).dx, 1) )
      disp("ERROR: mismatched frame length in calc_d2")
      return
   end 

   dx2_mat       = mask.*( trackmat_xyl(ti).dx .* trackmat_xyl(ti).dx );
   dy2_mat       = mask.*( trackmat_xyl(ti).dy .* trackmat_xyl(ti).dy );

   dx_list = [ dx_list , dx2_mat(mask == 1) ];
   dy_list = [ dy_list , dy2_mat(mask == 1) ];
   
end

d2 = mean( (dx_list+dy_list)/2  )
% Done this way so that if there is a mismatched dimension between
% x and y, the error will be flagged.

