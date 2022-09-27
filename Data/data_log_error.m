function [err, maximum] = data_log_error(Orig, Approx, delta)

% Usage : [rms_error, max_error] = data_log_error(Orig, Approx, delta=1E-6)
%
% Computes the rms y max relative error given the original and
% the approximation for the function. The order IS important.
% The value of delta is added to each entry to avoid log domain errors.

if nargin==2
  delta = 1E-6;
end

Orig = log(delta + reshape(Orig, [], 1));
Approx = log(delta + reshape(Approx, [], 1));
Diff = abs(Orig - Approx);
err = norm(Diff) / sqrt(size(Diff,1));
maximum = max(Diff);

end
