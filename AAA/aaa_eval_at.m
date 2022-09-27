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

% usage: h = aaa_eval_at(ss, fs, ws, freqs)
% Input: 
%   ss, fs, ws = vectors of support pts, function values, weights
%   freqs = Frequencies where to evaluate.
% Output:
%   h = function values at points freqs.

function hs = aaa_eval_at(freqs, ss, fs, ws)

sks = 2j * pi * freqs(:);
CC = 1 ./ bsxfun(@minus, sks, ss.');
hs = (CC * (ws.*fs)) ./ (CC * ws);
hs = reshape(hs, size(sks));

end
