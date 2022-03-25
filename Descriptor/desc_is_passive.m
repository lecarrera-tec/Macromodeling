function answer = desc_is_spassive(SS)

  % Usage : answer = desc_is_spassive(SS)
  %
  % Check if SS scatering parameters satisfies
  % passivity.
  
  assert(norm(SS.D) < 1);
  n = size(SS.A,1);
  II = eye(size(SS.D,1));
  [L, U, P] = lu(II - SS.D' * SS.D);
  assert(n == size(SS.A, 2));
  M = zeros(2 * n, 2 * n);
  M(1:n,1:n) = SS.A + SS.B * (U \ (L \ (P * SS.D'))) * SS.C;
  T = U \ (L \ (P * SS.B'));
  M(1:n,n+1:end) = SS.B * T;
  M(n+1:end, n+1:end) = -SS.A' - SS.C' * SS.D * T;
  M(n+1:end,1:n) = -SS.C' * ((II - SS.D * SS.D') \ SS.C);

  lambda = real(eig(M));
  answer = nnz(lambda) == size(lambda);
end
