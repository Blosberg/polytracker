function [ death_evi ] = find_death_event( SoE_in, tracknum )
% get the event index for the birth of tracknum

death_evi = find( SoE_in(:,2) == 2 && SoE(:,3) == tracknum );

if ( length(birth_evi)~=1 )
    disp("ERROR: non-unique death event");
end

end
