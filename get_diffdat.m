function [ Diffdat ] =  get_diff_dat ( statemat, trackmat_xyl, max_state, dt, Nframes )
% Take x-y track data,  

Num_comp_tracks = length( tracks_input);

for S= 1:max_state
    Diffdat{S} = calculate_diff_const(statemat, trackmat_xyl, Nframes, S);

    dx2                = calculate_dx2(statemat, trackmat_xyl, S);
    dxndxnp1           = calculate_dxnnp1(statemat, trackmat_xyl, S);

    % this is the diffusion constant, using the formula from Eq. 14 (page 022726-7
    % from Vestergaard et al, Phys. Rev. E. 89, (2014)
    Diffdat{S}.D       = (dx2/2*(dt)) + (dxndxnp1/dt)/dt;
 
%   R = 1/6 (consult with Jan about this: 
%   Diffdat{S}.sigma2  = R*(dx2) + (2*R-1)*dxndxnp1;
end


end

