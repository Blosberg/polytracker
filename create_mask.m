function [ mask_out ] = get_state_mask( state_matrix, desired_state, boolderiv )
% takes a matrix of states (state_matrix) and outputs a mask of 1's and 0's
% the ones indicate frames where the state_matrix == desired_state 
% if deriv ==0, then mask_out has the same size as state_matrix
% if deriv ==1, then one frame is removed, and the mask value == 1 only iff
% the state persists through the current from *and* the frame immediately 
% thereafter.
 
@@@

end

