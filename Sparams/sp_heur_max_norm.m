function [omega_max, hinf] = sp_heur_max_norm(omegas, Sys)
  % Usage: [omega_max, hinf] = sp_heur_max_norm(omegas, Sys)
  % A heuristic that look for the maximum norm in the
  % State Space or Descriptor form system representation.
  %
  % omegas : 2*pi*f, where f are a set of frequencies
  %          for the frequency range of interest.
  % Sys : State Space or descriptor form system representation.
  %
  % omega_max : Omega for maximum obtained norm.
  % hinf : Maximum obtained norm.

  % Computing the norm for each frequency.
  m = size(omegas, 2);
  norms = zeros(1,m);
  Hs = desc_eval(omegas, Sys);
  [nrows, ncols] = size(Sys.D);
  for ii = 1:m
    Hi = reshape(Hs(:,:,ii), nrows, ncols);
    norms(ii) = norm(Hi);
  end

  omega_max = 0;
  hinf = 0;
  for ii = 3:m
    xs = omegas(ii-2:ii);
    ys = norms(ii-2:ii);
    if max(ys) < 0.95 * hinf
      continue;
    end
    [xi, yi] = find_max(xs, ys, Sys);
    if hinf < yi
      omega_max = xi;
      hinf = yi;
    end
  end
end

