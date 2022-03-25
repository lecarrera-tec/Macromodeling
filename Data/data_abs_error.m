function [err, maximum] = data_abs_error(Orig, Approx)

  % Usage : [rms_error, max_error] = data_abs_error(Orig, Approx)
  %
  % Computes the absolute error given the original and
  % the approximation for the function.
  % The order IS important.
  
  Orig = reshape(Orig, [], 1);
  Approx = reshape(Approx, [], 1);
  Diff = abs(Orig - Approx);
  err = norm(Diff) / sqrt(size(Diff,1));
  maximum = max(Diff);
end
