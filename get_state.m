function [ state_output, window_duration ] =  get_state ( state_vector_ti_sti, Event_i, Event_f )
% having identified a process of interest, obtain the state at the time, and the duration.

% ---- sanity checks:
if( size( Event_i, 2 ) ~= 4 || size( Event_f, 2 ) ~= 4 )
  disp("ERROR: incorrectly sized Event inputs")
  return
end

if( Event_f(1) - Event_i(1) <=0)
  disp("ERROR: time window of process non-finite")
  return
end

temp_state_vector = state_vector_ti_sti( Event_i(1):Event_f(1)-1 );
S=unique( temp_state_vector( ~isnan(temp_state_vector) ) );

if( size(S,2) ~= 1 )
   disp("ERROR: non-uniform states during excitation window")
   return
end

if( state_vector_ti_sti( Event_i(1) ) -  state_vector_ti_sti( Event_i(1)-1 ) <=0  || state_vector_ti_sti( Event_f(1) ) - state_vector_ti_sti( Event_f(1)-1 ) >=0 )
   disp("ERROR: unexpected state transition trajectory during excitation window")
   return
end

state_output          = S;
% @@@ TODO: double-check +1 convention here:
window_duration = Event_f(1) - Event_i(1);

end
