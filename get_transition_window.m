function [ lifetime_list ] = get_state_lifetimes( state_matrix, Seq_of_events )
%  returns an array of all successive lifetimes (in units of frames)
%  for polymers of this state, using the sequence of events --
%  this is necessary because sometimes dimers/trimers split and then 
%  die in the same frame, so you can't just use changes in the state matrix


end

