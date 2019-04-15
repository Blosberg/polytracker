function [ fullDim_array ] =  convert_diffstruct2array( Struct_array_in )
% convert struct elements into single array --only valid for a very
% specific use-case.

N             = size( Struct_array_in, 1);
fullDim_array = [];

for i = 1:N
    fullDim_array = [ fullDim_array; Struct_array_in(i).fullDim.normDiffCoef( ~isnan( Struct_array_in(i).fullDim.normDiffCoef ) ) ];
end

end

