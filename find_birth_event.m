function [ birth_evi ] = find_birth_event( SoE_in, tracknum )
% get the event index for the birth of tracknum

% two sets of booleans get evaluated to 1/0, and then .* multiplied for AND
birth_evi = find( (SoE_in(:,2) == 1 ).*( SoE_in(:,3) == tracknum )  );

if ( length(birth_evi)~=1 )
    disp("ERROR: non-unique birth event");
end

end
