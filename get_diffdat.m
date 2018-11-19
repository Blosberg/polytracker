function [ Diffdat, dx2, dxndxnp1 ] =  get_diffdat ( statemat, trackdat_xyl, max_state, dt, Nframes ,R )
% Take x-y track data,

Num_comp_tracks = length( trackdat_xyl);

for S= 1:max_state

    dx2               = calculate_dx2(    statemat, trackdat_xyl, Nframes, S);
    dx2_ave           = mean(dx2);

    dxndxnp1          = calculate_dxnnp1( statemat, trackdat_xyl, Nframes, S);
    dxndxnp1_ave      = mean(dxndxnp1);

    % this is the diffusion constant, using the formula from Eq. 14 (page 022726-7
    % from Vestergaard et al, Phys. Rev. E. 89, (2014)
    Diffdat{S}.D       = (dx2_ave/2*(dt)) + (dxndxnp1_ave/dt)/dt;
    Diffdat{S}.sigma2  = R*(dx2_ave) + (2*R-1)*dxndxnp1_ave;

    if( Diffdat{S}.sigma2 <1  )
        disp( strcat("WARNING: sigma^2 =", num2str(Diffdat{S}.sigma2 ), " for state ",num2str( S ) ) )
    end

end

%-------- now get the overall (stateless) diff constant distribution
%-------- (i.e. calculated for each subtrack)

% %%%%%% START HERE @@@@@@
%
% for ti = 1:Num_comp_tracks
%
%    % get mask for frames where the track was in the right state in that
%    % frame AND the subsequent one:
%    mask          =  create_mask( statemat, Nframes, 0, 1 );
%
%    if( ~all( size(mask) == [size( trackdat_xyl(ti).dx, 1), (Nframes-1) ] ) )
%       disp("ERROR: mismatched frame length in calc_d2")
%       return
%    end
%
%    dx2_mat       = mask.*( trackdat_xyl(ti).dx .* trackdat_xyl(ti).dx );
%    dy2_mat       = mask.*( trackdat_xyl(ti).dy .* trackdat_xyl(ti).dy );
%
%    mean(dx2_mat + dy2_mat)'
%
%    dx2_vec = dx2_mat(mask == 1);
%    dy2_vec = dy2_mat(mask == 1);
%
%    % if there's only a single subtrack, then the array gets returned with a different orientation.
%    % Make sure it's a horizonal array:
%    if( size(dx2_vec,1) ~=1  )
%        dx2_vec = transpose( dx2_vec );
%        dy2_vec = transpose( dy2_vec );
%    end
%
%    dx_list = [ dx_list , dx2_vec ];
%    dy_list = [ dy_list , dy2_vec ];
%
% end
% %%%%%%
%
end

