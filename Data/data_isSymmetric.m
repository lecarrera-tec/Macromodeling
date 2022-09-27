function p = data_isSymmetric(Hs, tol)
% Usage: p = data_isSymmetric(Hs, tol = 1e-6)
% Compute for each transfer function if it
% is symmetric, using a relative tolerance tol

if nargin == 1
  tol = 1e-6;
end

p = false;
nports = size(Hs, 1);
temp = size(Hs, 2);
if nports ~= temp
  return;
end

nfreqs = size(Hs, 3);
for ii = 1:nfreqs
  Hi = reshape(Hs(:,:,ii), nports, nports);
  temp = norm(Hi - transpose(Hi), 'fro') / norm(Hi, 'fro');
  if temp > tol
    return;
  end
end

p = true;
end
