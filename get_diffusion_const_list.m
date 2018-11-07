function [ lumen_list ] =  get_diffusion_const_list ( tracks_input, state_matrix_allti, tracksmat_xyl, max_state , Nframes )
% builds a list of observed luminescance intensities for each type of polymer (i.e. size{1}, {2}, etc.)

Num_comp_tracks     = length( tracks_input );

% initialize the list
for S= 1:max_state
    lumen_list{S} = [];
end
% ---------------------------------

for ti = 1:Num_comp_tracks

  Nsubtracks(ti) =  size( tracks_input(ti).tracksCoordAmpCG, 1);
     
  for S = 1:max(max(state_matrix_allti{ti}))   
     
      mask             =  create_mask( state_matrix{ti}, S, 0 ) 
        
      dx_temp       = resize( mask.*dx{ti} ) @@@ grab dx data above
      Polydat{S}.dx = [ Polydat{S}.Lumen, Lumen_temp[ maskA != 0 ] ]

      dy_temp       = resize( mask.*dy{ti} ) @@@ grab dy data above
      Polydat{S}.dy = [ Polydat{S}.Lumen, Lumen_temp[ maskA != 0 ] ]
      
   end
end