function [ comptrack_lengths ] =  get_comptrack_lengths( tracksFinal_in )
% obtain an array of the time-lengths of composite tracks.

Ncomptracks  = size(tracksFinal_in);

for ti = 1:Ncomptracks
   Nentries =  size(tracksFinal_in(ti).tracksCoordAmpCG, 2);
   if ( mod (Nentries,8) ~= 0)
       print("ERROR: number of frames isn't a multiple of 8, somethings wrong here!!")
   end

   comptrack_lengths(ti) = size(tracksFinal_in(ti).tracksCoordAmpCG,2)/8;
end

comptrack_lengths = comptrack_lengths';


end
