function [freqs, Hs] = data_read(filename)

% Usage : [freqs, Hs] = data_read(filename)
%
% Reads a .data file of a square system of >nports<
% size, and extracts the frequency vector
% of size >nfreqs< and the transfer functions as a 3D matrix
% of size (nports, nports, nfreqs).
%
% The first line of the data should be:
% <nfreqs> <nports> <scale>
%
% Returns the vector of frequencies 
% and the transfer functions.

fid = fopen(filename, 'r');
infile = fscanf(fid, '%lf %lf %lf', [1,3]);
nfreqs = infile(1);
nports = infile(2);
scale  = infile(3);
nfuncts = nports * nports;
freqs = zeros(1, nfreqs);
Hs = zeros(nports, nports, nfreqs);
for k = 1:nfreqs
  freqs(k) = fscanf(fid, '%lf', [1,1]);
  row  = fscanf(fid, '%lf %lf', [2,nfuncts]);
  row = row';
  temp = reshape(row(:,1) + 1i*row(:,2), nports, nports);
  Hs(:,:,k) = temp.';
end
fclose(fid);

if (scale ~= 1)
  freqs = scale * freqs;
end

end
