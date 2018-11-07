function [ mask_out ] = create_mask( state_matrix, desired_state, boolderiv )
% takes a matrix of states (state_matrix) and outputs a mask of 1's and 0's
% the ones indicate frames where the state_matrix == desired_state 
% if deriv ==0, then mask_out has the same size as state_matrix
% if deriv ==1, then one frame is removed, and the mask value == 1 only iff
% the state persists through the current from *and* the frame immediately 
% thereafter.

if ( ~boolderiv )
    % -- considering only current state
    mask_out = state_matrix;
    mask_out( mask_out ~= desired_state ) = 0;
    mask_out( mask_out == desired_state ) = 1;
elseif ( boolderiv )

    init_state_match = state_matrix(:,1:(size(state_matrix,2)-1) ) ;
    % trim the last state to get the state at the beginning of the window 
    init_state_match( init_state_match ~= desired_state ) = 0;
    init_state_match( init_state_match == desired_state ) = 1;

    %  Y = diff(X,n,dim) is the nth difference calculated along the dimension specified by dim. The dim input is a positive integer scalar.    
    temp     = diff( state_matrix, 1, 2); 
    mask_out = ( (~temp) .* init_state_match  );
    % The particle was in the desired state at the onset of the
    % window, and had zero change at the end of the window.
end
