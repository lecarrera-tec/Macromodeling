function [err, maximum] = data_log_error(Orig, Approx, delta=1e-6)

  % Usage : [rms_error, max_error] = data_log_error(Orig, Approx)
  %
  % Computes the rms y max relative error given the original and
  % the approximation for the function.
  % The order IS important.

  Orig = log(delta + reshape(Orig, [], 1));
  Approx = log(delta + reshape(Approx, [], 1));
  Diff = abs(Orig - Approx);
  err = norm(Diff) / sqrt(size(Diff,1));
  maximum = max(Diff);
end
