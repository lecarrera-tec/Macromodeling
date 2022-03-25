function [DF, n] = desc_stab_flip(DF)

  % Usage : DF = desc_stab_flip(DF)
  % Given a Descriptor System DF, after transforming
  % the System to Upper Triangular matrices A and E,
  % it flips the sign of the real part for the unstable
  % poles.

  global fid;

  [A, E, Q, Z] = qz(DF.A, DF.E);
  Lambdas = ordeig(A, E);
  DF.B = Q * DF.B;
  DF.C = DF.C * Z;

  Idx = find(real(Lambdas) >= 0);
  Lambdas = Lambdas(Idx);
  n = length(Idx);
  fprintf(fid, "Number of unstable poles = %d\n", n);

  for ii = 1:n
    if imag(Lambdas(ii)) == 0
      E(Idx(ii), Idx(ii)) = -E(Idx(ii), Idx(ii));
    else
      A(Idx(ii), Idx(ii)) = -A(Idx(ii), Idx(ii));
    end
  end

  DF.A = A;
  DF.E = E;
end
