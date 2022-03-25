function [omegas, sigmas]  = sp_passv(DF, omax, opts)

% Usage: [omegas, sigmas] = sp_passv(DF, omax, opts)
%
% Gives omegas where there exists passivity violations.
% In the list is included hinf. Look for the greatest sigma.
% State Space or Descriptor form system representation.
%
% DF     : State Space or descriptor form system representation.
% omax    : Max frequency to look for. Starts from 0.
% opts.tol      : Tolerance for cero (default is 1). If <= 0, 
%                 uses only random points.
% opts.nsamples : Number of total samples (500 by default).
% opts.first    : Is the first time with random samples. Use linspace.
%               : Default is false.
%
% omegas : Omegas.
% sigmas : Respective values of sigmas greater that one.
%
% The steps of the algorithm are:
% 1. Compute passivity violations L using Hamiltonian matrix.
%    (all the elements of L have norm=1).
% 2. Make a list S = linspace(0, omax, opts.nsamples)
% 3. Compute the norms in the M points: at the beggining in
%    omega = 0, midpoint of L and omega = omax.
% 4. If between L(i) and L(i+1) is not passive, add the elements
%    of S that belong to such interval to a list K. Compute the
%    norm of each element of K.

if nargin == 2
  opts = 0;
end

if isfield(opts, 'nsamples')
  nsamples = opts.nsamples;
else
  nsamples = 500;
end

if isfield(opts, 'tol')
  tol = opts.tol;
else
  tol = 1;
end

if isfield(opts, 'first')
  firsttime = opts.first;
else
  firsttime = false;
end

[nrows, ncols] = size(DF);

% 1. Computing passivity violations by Hamiltonian matrix.
if tol > 0
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

  L = imag(L(abs(real(L)) < tol));
  L = L(L >= 0);
  L = sort(L(L <= omax));
  if omax == Inf
    omax = 1.1 * L(end);
  end
  clear MM K0

  if isempty(L)
    omegas = [];
    sigmas = [];
    return
  end

  % 2. Full list of omega samples.
  omegas = L;
  sigmas = ones(size(L));
  S = linspace(0, omax, nsamples);

  % 3. Extreme and midpoints where to do the first test
  M = [0, 0.5 * ( L(1:end-1) + L(2:end) ), omax];
  m = size(M, 2);
  norms = desc_norms(DF, M);
  omegas = [omegas, M];
  sigmas = [sigmas, norms];

  % Non-passive at omega = 0.
  if sigmas(1) >= 1
    list = S(logical((S > 0) .* (S < L(1))));
  else
    list = [];
  end

  % Checking if the model is passive or not between 
  % the points given by the Hamiltonian matrix.
  for ii = 2:m-1
    if sigmas(ii) >= 1 
      list = [list, S(logical((S > L(ii - 1)) .* (S < L(ii))))];
    end
  end

  % Checking passivity violations at the end.
  if sigmas(end) >= 1
    list = [list, S(logical((S > L(end)) .* (S < omax)))];
  end
  norms = desc_norms(DF, list);

  % We add and sort the values already computed.
  [omegas, idx] = sort([omegas, list]);
  sigmas = [sigmas, norms];
  sigmas = sigmas(idx);
elseif firsttime % We are not using a Hamiltonian matrix. 
    % Use a linear sample points.
    omegas = linspace(0, omax, nsamples);
  sigmas = desc_norms(DF, omegas);
else % Use a random matrix instead.
  omegas = sort([0, omax * rand(1, nsamples), omax]);
  sigmas = desc_norms(DF, omegas);
end % if tol > 0

% Looking for hinf using a quadratic vertex.
oldmax = 0;
smax = max(sigmas);
while smax > oldmax
  oldmax = smax;

  list = [];
  m = length(omegas);
  for ii = 2:m-1
    if sigmas(ii) > sigmas(ii-1) && sigmas(ii) > sigmas(ii+1)
      LHS = [omegas(ii-1); omegas(ii); omegas(ii+1)] .^ [2, 1, 0];
      RHS = [sigmas(ii-1); sigmas(ii); sigmas(ii+1)];
      X = LHS \ RHS;
      list(end+1) = -0.5 * X(2) / X(1);
    end
  end
  norms = desc_norms(DF, list);
  [omegas, idx] = sort([omegas, list]);
  sigmas = [sigmas, norms];
  sigmas = sigmas(idx);
  [omegas, sigmas] = del_repeated(omegas, sigmas);
  smax = max(sigmas);
end

idx = sigmas >= 1;
sigmas = sigmas(idx);
omegas = omegas(idx);

end
