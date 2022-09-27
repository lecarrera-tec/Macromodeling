function Hs = desc_eval(freqs, Sys)
% Usage : Hs = desc_eval(freqs, Sys)
% 
% Computes the transfer function, given the
% (generalized) space state representation 
% and desired frequency points as \omega.

jomegas = 2j * pi * freqs;

if isfield(Sys, 'E')
  [A, E, Q, Z] = qz(Sys.A, Sys.E);
  B = Q * Sys.B;
  C = Sys.C * Z;
else
  A = Sys.A;
  E = eye(size(A));
  B = Sys.B;
  C = Sys.C
end; %if isfield(Sys, 'E')

nrows = size(C, 1);
ncols = size(B, 2);
nfreqs = length(omegas);
Hs = zeros(nrows, ncols, nfreqs);

for k = 1:nfreqs
  % Solve the system.
  X = jomegas(k) * E - A;
  Hs(:,:,k) = C * (X \ B) + Sys.D;
end;

end
