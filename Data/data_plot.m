function data_plot(freqs, Hs, filename)
% Usage data_plot(freqs, Hs, filename)
%
% Saves in the filename the magnitud in db 
% and the argument in radians. Hs is a transfer
% function for the given frequencies.

isOctave = (exist('OCTAVE_VERSION') > 0);

nfreqs = size(freqs,2);
Mags  = 20 * log10(abs(Hs));
Args  = angle(Hs);
Nrows = size(Hs,1);
Ncols = size(Hs,2);

% We need to include the frequencies, the magnitudes
% and the arguments.
P = zeros(nfreqs, 1 + 2 * Nrows * Ncols);

P(:,1) = freqs';
index = 2;
for irow = 1:Nrows
  for icol = 1:Ncols
    P(:, index) = reshape(Mags(irow, icol, :), [], 1);
    index = index + 1;
    P(:, index) = reshape(Args(irow, icol, :), [], 1);
    index = index + 1;
  end % for icol = 1:Ncols
end % for irow = 1:Nrows

if isOctave
  save('-ascii', filename, 'P');
else
  save(filename, 'P', '-ascii');
end

end % endfunction
