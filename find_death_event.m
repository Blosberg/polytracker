function [ death_evi ] = find_death_event( SoE_in, tracknum )
% get the event index for the birth of tracknum

% two sets of booleans get evaluated to 1/0, and then .* multiplied for AND
death_evi = find( (SoE_in(:,2) == 2 ).*( SoE_in(:,3) == tracknum ) );

if ( length(death_evi)~=1 )
    disp("ERROR: non-unique death event");
end

end
