% Copyright (C) Nataksukasa, Sete, Trefethen
%
% This is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 3 of the License, or (at your option) any
% later version.
%
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details: <http://www.gnu.org/licenses/>.
%
% Adaptaded-By: Luis Ernesto Carrera-Retana <lecarrera@itcr.ac.cr>
% Following and idea from FastAAA: A fast rational-function
% fitter by Amit Hochman.

% aaa rational approximation of data hs on set ss
% [s, f, w, ds] = aaa(freqs, hs, tol, mmax, rel_error)
%
% Input: freqs = vector of frequency sample points
%        hs = vector of data values, or a function handle
%        tol = tolerance tol, set to 1e-13 if omitted
%        mmax: max type is (mmax-1,mmax-1), set to 100 if omitted
%        rel_error: use relative instead of absolute error.
%
% Output: 
%     s,f,w = vectors of support pts, function values, weights
%     ds    = vector of ds values on each iteration.
%     ers   = vector of errors on each iteration.

function [s, f, w, ds, ers] = aaa_realfit(freqs, hs, tol, mmax, rel_error)

global fid;

ss = 2j * pi * freqs;

% default tol 1e-13
if nargin<3, tol = 1e-13; end
% default max type (99,99)
if nargin<4, mmax = 100; end
% default error
if nargin<5, rel_error = false; end

ss = vec(ss, 2);
ss = vec([ss; -ss]);
hs = vec(hs, 2);
hs = vec([hs; conj(hs)]);

% number of sample points
M = length(ss);

% left scaling matrix
SF = spdiags(hs,0,M,M);

% initializations
J = 1:M; s = []; f = []; C = [];
R = mean(hs);
ds = [];

if rel_error
  stop_criterion = tol;
else
  stop_criterion = tol*norm(hs, inf);
endif;

% main loop
for m = 1:mmax
  % select next support points
  if rel_error
    [~, jj] = max(abs(hs-R) ./ abs(hs));
  else
    [~,jj] = max(abs(hs-R));
  endif;
  if jj > 1 && ss(jj-1) == -ss(jj)
    jj -= 1;
  endif;
  % update support points, data values
  s = [s; ss(jj); ss(jj+1)]; f = [f; hs(jj); hs(jj+1)];
  % update index vector
  J(J==jj) = [];
  J(J==jj+1) = [];
  % next two columns of Cauchy matrix
  C = [C, 1./(ss-ss(jj)), 1 ./ (ss - ss(jj+1))];
  % right scaling matrix
  Sf = diag(f);
  % Loewner matrix
  A = SF*C - C*Sf;

  %%%%%%%%%%%%%%%%%%%%%%
  % Hermitian symmetry %
  %%%%%%%%%%%%%%%%%%%%%%
  LH = A(J(1:2:end),:);
  [nrows, ncols] = size(LH);
  LH = reshape([real(LH); imag(LH); -imag(LH); real(LH)], 2 * nrows, 2 * ncols);
  H = kron(eye(0.5 * ncols), [1, 0; 0, 1; 1, 0; 0, -1]);

  % SVD
  [~,~,V] = svd(LH * H,0);

  % weight vector = min sing vector
  w = H * V(:, end);
  w = w(1:2:end) + 1j * w(2:2:end);
  %%%%%%%%%%%%%%%%%%%%%%
  % Hermitian symmetry %
  %%%%%%%%%%%%%%%%%%%%%%

  % numerator and denominator
  N = C*(w.*f); D = C*w;
  % rational approximation
  R = hs; R(J) = N(J)./D(J);
  % max error at sample points
  if rel_error
    err = norm(abs(hs-R) ./ abs(hs), inf);
  else
    err = norm(hs-R,inf);
  endif;
  ds = [ds, vec(f, 2) * vec(w, 1) / sum(w)];
  % stop if converged
  if err <= stop_criterion, break, end
end % for m = 1:mmax

if err > stop_criterion
  error('real-AAA did not converge');
end

end % function
