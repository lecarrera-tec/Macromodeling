function ns = desc_norms(DF, omegas)
% Usage: ns = desc_norms(DF, omegas)
% Computes the maximum singular value
% of the descriptor system at each given omega.
% DF : Descriptor system.
% omegas : frequencies to compute the norms.

% Computing the norm for each new frequency.
m = size(omegas, 2);
ns = zeros(1, m);
[nrows, ncols] = desc_size(DF);
kw = 0;
jump = 21;
while kw < m
  k0 = kw;
  kw = min(m, k0 + jump);
  Hs = desc_eval(omegas(k0+1:kw), DF);
  for ii = k0+1:kw
    ns(ii) = norm(reshape(Hs(:,:,ii-k0), nrows, ncols));
  end
end

end
