% Copyright (C) <lecarrera at itcr dot ac dot cr>
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
% Author: Luis Ernesto Carrera-Retana

% Usage: [Hps, D] = aaa_full(freqs, Hs, opts)
% Computes the derivative of a reciprocal system Hs,
% and an approximation to the matrix D.
% Input:
%   opts.real  : Use a real fitting (default).
%   opts.dext : d extraction.
%     false : not extract.
%   opts.dreal     : if using a complex fitting, use the real part (default).
%                    If not, use the absolute value with the sign of the real part.
%   opts.rel_error : Use relative error (default). If not, use abs error.
%   opts.tol       : Tolerance used in the fitting algorithm.
%   opts.mmax      : Maximum number of iterations.
%   opts.horder    : High order derivative.
%   opts.ssize     : Step size for the derivative.

function [Hps, D] = aaa_full(freqs, Hs, opts)

global fid;

%-%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
%-%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
  opts.real = true;
endif;

if ~isfield(opts, 'real')
  opts.real = true;
endif;

if ~isfield(opts, 'dtype')
  opts.dtype = 'median';
endif;

if ~opts.real && ~isfield(opts, 'dreal')
  opts.dreal = true;
endif;

if ~isfield(opts, 'rel_error')
  opts.rel_error = true;
endif;

if ~isfield(opts, 'tol')
  opts.tol = 1e-6;
endif;

if ~isfield(opts, 'horder')
  opts.horder = true;
endif;

if ~isfield(opts, 'ssize')
  opts.ssize = 0.1;
endif;

nfreqs = length(freqs);
if ~isfield(opts, 'mmax')
  opts.mmax = floor(0.9 * nfreqs);
endif;

% Choosing algorithms
if opts.real
  fprintf(fid, 'Using real fitting ...\n');
  aaa = @aaa_realfit;
  dalg = @real_median_d;
else
  fprintf(fid, 'Using complex fitting ...\n');
  aaa = @aaa_fitting;
  dalg = false;
endif;

if opts.horder
  diff = @ihorder;
else
  diff = @ilorder;
endif;

nports = size(Hs, 1);
assert(nports == size(Hs, 2));
assert(nfreqs == size(Hs, 3));

% Initialization
Hps = zeros(size(Hs));
if ~isbool(dalg)
  D = zeros(nports);
endif;

for jj = 1:nports
  for ii = 1:jj-1
    hs = vec(Hs(ii, jj, :), 2);
    [s, f, w, ds] = aaa(freqs, hs, opts.tol, opts.mmax, opts.rel_error);
    Hps(ii,jj,:) = diff(s, f, w, freqs, opts.ssize);
    Hps(jj,ii,:) = Hps(ii,jj,:);
    if ~isbool(dalg)
      D(ii, jj) = dalg(ds, hs);
      D(jj, ii) = D(ii, jj);
    endif;
  endfor;
  hs = vec(Hs(jj, jj, :), 2);
  [s, f, w, ds] = aaa(freqs, hs, opts.tol, opts.mmax, opts.rel_error);
  Hps(jj,jj,:) = diff(s, f, w, freqs, opts.ssize);
  if ~isbool(dalg)
    D(jj, jj) = dalg(ds, hs);
  endif
endfor;

% Passivity
if ~isbool(dalg) & norm(D) >= 1
  [U, S, V] = svd(D);
  allowed = (1 - 1e-6) * ones(size(S));
  S = min(allowed, S);
  D = U * S * V';
endif;

endfunction;

function yp = horder(s, f, w, freqs, h)
  sks = 2j * pi * freqs;
  yp = aaa_eval_at(s, f, w, sks - 2i*h);
  yp -= 8 * aaa_eval_at(s, f, w, sks - 1i*h);
  yp += 8 * aaa_eval_at(s, f, w, sks + 1i*h);
  yp -= aaa_eval_at(s, f, w, sks + 2i*h);
  yp /= 12i * h;
endfunction;

function yp = lorder(s, f, w, freqs, h)
  sks = 2j * pi * freqs;
  yp = aaa_eval_at(s, f, w, sks + 1i*h);
  yp -= aaa_eval_at(s, f, w, sks - 1i*h);
  yp /= 2i * h;
endfunction;

function yp = ihorder(s, f, w, freqs, h)
  sks = 2j * pi * freqs;
  yp = aaa_eval_at(s, f, w, sks - 2*h);
  yp -= 8 * aaa_eval_at(s, f, w, sks - h);
  yp += 8 * aaa_eval_at(s, f, w, sks + h);
  yp -= aaa_eval_at(s, f, w, sks + 2*h);
  yp /= 12 * h;
endfunction;

function yp = ilorder(s, f, w, freqs, h)
  sks = 2j * pi * freqs;
  yp = aaa_eval_at(s, f, w, sks + h);
  yp -= aaa_eval_at(s, f, w, sks - h);
  yp /= 2 * h;
endfunction;

function d = real_median_d(ds, hs)
  dmin = max(-1, min(hs));
  dmax = min(1, max(hs));
  drng = max(-dmin, dmax);
  ds(ds < -drng) = [];
  ds(ds > drng) = [];
  d = median(ds);
endfunction;
