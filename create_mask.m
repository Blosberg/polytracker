function [ mask_out ] = create_mask( state_matrix, Nframes, desired_state, deriv )
% takes a matrix of states (state_matrix) and outputs a mask of 1's and 0's
% the ones indicate frames where the state_matrix == desired_state
% if deriv ==0, then mask_out has the same size as state_matrix
% if deriv ==1/2, then one frame is removed, and the mask value == 1 only iff
% the state persists through the current from *and* (1) the frame immediately
% thereafter, or (2) the next two frames.

if( size(state_matrix,2) ~= Nframes )
   disp("ERROR: missized state matrix")
   return
end


% if desired_state == 0, this is code for any finite state >0 (since there is no
% such state 0 for a polymer
if (desired_state ==0)
   % set every state >1 to 1, and set desired state to 1. Then produce the
   % monomer mask as usual
   state_matrix = (state_matrix >=1)
   desired_state =1
end
% ---


if ( deriv ==0 )
    % -- states must match in CURRENT frame
    mask_out          = state_matrix;

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
    state_change                       = diff( state_matrix, 1, 2);
    state_change( isnan(state_change)) = 1;
    mask_out           = ( (~state_change) .* init_state_match  );
    % The particle was in the desired state at the onset of the
    % window, and experienced no change in the following frame.

elseif ( deriv == 2 )
    % states must match in current frame AND the next TWO frames

    init_state_match = state_matrix(:,1:(Nframes-2) ) ;
    % trim the last 2 states to get the state at the beginning of the window

    init_state_match( init_state_match ~= desired_state ) = 0;
    % set any non-matching frames to 0

    init_state_match( init_state_match == desired_state ) = 1;
    % set any matching frames to 1

    %  Y = diff(X,n,dim) is the nth difference calculated along the dimension specified by dim. The dim input is a positive integer scalar.
    state_change1     = diff( state_matrix, 1, 2);
    state_change1( isnan(state_change1)) = 1;

    state_preserved2 = ~state_change1(:,1:(Nframes-2) ) .* ~state_change1(:,2:(Nframes-1) );
    % check whether the "current" state is preserved across the next *TWO* frames

    mask_out = ( (state_preserved2) .* init_state_match  );
    % The particle was in the desired state at the onset of the
    % window, and had zero change over the following *TWO* frames.
else
    disp("ERROR: invalid value submitted for deriv ")
    return
end
