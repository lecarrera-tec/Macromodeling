function m = desc_pass_measure(DF, omax, use_hamil=true)

% Usage: m = desc_pass_measure(DF, omax, use_hamil = true)
% Computes the area of the function \sum max(0, \sigma_i-1)
% for every singular value.

global fid;
fprintf(fid, 'Measuring passivity violations:');

if nargin == 2
  use_hamil = true
end

k = 0.5e-9 / pi;
[nrows, ncols] = desc_size(DF);
tt = omax * nrows;
m = [0, 0];

f = @(X) one_sigma(DF, X);
g = @(X) all_sigmas(DF, X);

if use_hamil
  SS = eye(nrows) - DF.D * DF.D';
  RR = eye(ncols) - DF.D' * DF.D;
  TT = SS \ DF.C;
  M11 = DF.A + DF.B * DF.D' * TT;
  M21 = -DF.C' * TT;
  TT = RR \ DF.B';
  M12 = DF.B * TT;
  M22 = -DF.A' - DF.C' * DF.D * TT;
  MM = [M11, M12; M21, M22];
  clear II RR M11 SS M12 M22 M21 M22

  K0 = blkdiag(DF.E, DF.E');
  L = transpose(eig(MM, K0));
  L = imag(L(abs(real(L)) < 1));
  L = sort(L(L > 0 & L < omax));

  clear SS RR TT M11 M21 TT M12 M22 MM K0
end

if ~use_hamil || isempty(L)
  ys = 0.01 * omax * (0:9).^1.7;
  ys = [ys, 0.5 * omax, flip(omax - ys)];
else
  ys = [0, L, omax];
end
  
nfreqs = length(ys) - 1;
for ii = 1:nfreqs
  pp = 100 * ii / nfreqs;
  temp = quadgk(f, ys(ii), ys(ii+1));
  if temp > 0
    m(1) = m(1) + temp / omax;
    m(2) = m(2) + quadgk(g, ys(ii), ys(ii+1)) / tt;
  end
  fprintf(fid, '; (%.2f-%.2f, %.2e/%.2e) %.1f%%', k * ys(ii), k * ys(ii+1), m(1), m(2), pp);
end
fprintf(fid, '\n');

end

function ys = all_sigmas(DF, X)
  global fid;
  ys = zeros(size(X));
  Hs = desc_eval(X, DF);
  n = size(Hs, 3);
  assert(n == length(X));
  [nrows, ncols] = desc_size(DF);
  for ii = 1:n
    Hi = reshape(Hs(:,:,ii), nrows, ncols);
    sigmas = svd(Hi);
    ys(ii) = sum(sigmas(sigmas > 1) - 1);
  end
end

function ys = one_sigma(DF, X)
  global fid;
  ys = zeros(size(X));
  Hs = desc_eval(X, DF);
  n = size(Hs, 3);
  assert(n == length(X));
  [nrows, ncols] = desc_size(DF);
  for ii = 1:n
    Hi = reshape(Hs(:,:,ii), nrows, ncols);
    ys(ii) = max(1, norm(Hi)) - 1;
  end
end
