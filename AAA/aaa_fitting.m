%% Copyright (C) Nataksukasa, Sete, Trefethen
%%
%% This is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public
%% License as published by the Free Software Foundation;
%% either version 3 of the License, or (at your option) any
%% later version.
%%
%% This code is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied
%% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
%% PURPOSE.  See the GNU General Public License for more
%% details: <http://www.gnu.org/licenses/>.
%%
%% Adaptaded-By: Luis Ernesto Carrera-Retana <lecarrera@itcr.ac.cr>

%% usage: [s,f,w] = aaa(hs,ss,tol,mmax,rel_error)
%% 
%% aaa rational approximation of data hs on set ss
%%
%% Input: hs = vector of data values, or a function handle
%%        ss = vector of sample points
%%        tol = tolerance tol, set to 1e-13 if omitted
%%        mmax: max type is (mmax-1,mmax-1), set to 100 if omitted
%%
%% Output: 
%%     s,f,w = vectors of support pts, function values, weights
%%     ds    = vector of ds values on each iteration.

function [s, f, w, ds] = aaa_fitting(freqs, hs, tol, mmax, rel_error)

global fid;

% default tol 1e-13
if nargin<3, tol = 1e-13; end
% default max type (99,99)
if nargin<4, mmax = 100; end
% default error
if nargin<5, rel_error = false; end

ss = vec(2j * pi * freqs);
hs = vec(hs);

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
  % update support points, data values
  s = [s; ss(jj)]; f = [f; hs(jj)];
  % update index vector
  J(J==jj) = [];
  % next column of Cauchy matrix
  C = [C, 1./(ss-ss(jj))];
  % right scaling matrix
  Sf = diag(f);
  % Loewner matrix
  A = SF*C - C*Sf;
  % SVD
  [~,~,V] = svd(A(J,:),0);
  % weight vector = min sing vector
  w = V(:, m);
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
  % stop if converged
  if err <= stop_criterion, break, end
end % for m = 1:mmax

if err > stop_criterion
  error('AAA did not converge');
end

ds = 0;

end % function
