function [ mask_out ] = create_mask( state_matrix, desired_state, deriv )
% takes a matrix of states (state_matrix) and outputs a mask of 1's and 0's
% the ones indicate frames where the state_matrix == desired_state 
% if deriv ==0, then mask_out has the same size as state_matrix
% if deriv ==1, then one frame is removed, and the mask value == 1 only iff
% the state persists through the current from *and* the frame immediately 
% thereafter.


if( size(state_matrix,2) ~= Nframes )
   disp("ERROR: missized state matrix")
   return
end

% ---


if ( deriv ==0 )
    % -- states must match in CURRENT frame
    mask_outi          = state_matrix;

    mask_out( mask_out ~= desired_state ) = 0;
    % set any non-matching frames to 0

    mask_out( mask_out == desired_state ) = 1;
    % set any matchine frames to 1

elseif ( deriv == 1 )
    % -- states must match in current AND subsequent frame

    init_state_match = state_matrix(:,1:(Nframes-1) ) ;
    % trim the last state to get the state at the beginning of the window 

    init_state_match( init_state_match ~= desired_state ) = 0;
    % set any non-matching frames to 0

    init_state_match( init_state_match == desired_state ) = 1;
    % set any matching frames to 1

    %  Y = diff(X,n,dim) is the nth difference calculated along the dimension specified by dim. The dim input is a positive integer scalar.    
    temp     = diff( state_matrix, 1, 2); 
    mask_out = ( (~temp) .* init_state_match  );
    % The particle was in the desired state at the onset of the
    % window, and had zero change at the end of the window.

elseif ( deriv == 2 )
    % states must match in current frame AND the next TWO frames

    init_state_match = state_matrix(:,1:(Nframes-2) ) ;
    % trim the last 2 states to get the state at the beginning of the window 

    init_state_match( init_state_match ~= desired_state ) = 0;
    % set any non-matching frames to 0

    init_state_match( init_state_match == desired_state ) = 1;
    % set any matching frames to 1

    %  Y = diff(X,n,dim) is the nth difference calculated along the dimension specified by dim. The dim input is a positive integer scalar.    
    temp     = ~diff( state_matrix, 1, 2);
    temp2    = temp( :, 1:(Nframes-2) ) .* temp( :, 2:(Nframes-1) )  
 
    mask_out = ( (temp) .* init_state_match  );
    % The particle was in the desired state at the onset of the
    % window, and had zero change at the end of the window.
else
    disp("ERROR: invalid value submitted for deriv ")
    return
end
