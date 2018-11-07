function [ lumen_list ] =  get_lumen_list ( tracks_input, state_matrix, tracksmat_xyl, max_state , Nframes )
% builds a list of observed luminescance intensities for each type of polymer (i.e. size{1}, {2}, etc.)

Num_comp_tracks     = length( tracks_input );

% initialize the list
for S= 1:max_state
    lumen_list{S} = [];
end
% ---------------------------------

for ti = 1:Num_comp_tracks

  max_state_for_this_comptrack = max(max(state_matrix{ti}));
   
  for S = 1:max_state_for_this_comptrack

     mask                  =  create_mask( state_matrix{ti}, S, 0 );
     if ( ~all( size(mask) == size( tracksmat_xyl(ti).Lamp ) ) )
         disp("mismatched mask dimensions for Lamp")
         return
      end    
    
    if (  size(mask,1) ==1)       
      lumen_list{S}    = [ lumen_list{S} ,  tracksmat_xyl(ti).Lamp(mask==1) ];    
    else
      lumen_list{S}    = [ lumen_list{S} ,  transpose(tracksmat_xyl(ti).Lamp(mask==1)) ]; 
    end
   end
   
end