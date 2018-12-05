function [ birth_evi ] = find_birth_event( SoE_in, tracknum)
% get the event index for the birth of tracknum

birth_evi = find( SoE_in(:,2) == 1 && SoE(:,3) == tracknum );

if ( length(birth_evi)~=1 )
    disp("ERROR: non-unique birth event");
end

end
