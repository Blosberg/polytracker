function [ tenable ] =  check_states_tenable ( SoE_statemat_in )
% returns 1 iff all numeric values in SoE_in are ge 1

tenable = min( SoE_statemat_in( ~isnan( SoE_statemat_in ) ) ) >= 1;

end
