function Sys = desc_updateB(Sys, Hs, omegas)
  % Usage: Sys = desc_updateB(Sys, Hs, omegas)
  % Update B to minimize the error.
  % It is an overdetermined system of the form
  % LHS * B = RHS, for the real and imaginary part.

  % Construct first matrices LHS and RHS.
  global fid;
  fprintf(fid, 'Updating B ...\n');

  nfreqs = size(omegas, 2);
  ssize = size(Sys.A, 1);
  [nrows, ncols] = desc_size(Sys);
  step = nrows * nfreqs;
  SysPhi = Sys;
  SysPhi.B = eye(ssize);
  SysPhi.D = 0;
  Phi = desc_eval(omegas, SysPhi);

  indi = 0;
  indf = 0;
  LHS = zeros(2 * step, ssize);
  RHS = zeros(2 * step, ncols);
  for i = 1:nfreqs
    indi = indf + 1;
    indf = indf + nrows;
    Temp = reshape(Phi(:, :, i), nrows, ssize);
    Hi = reshape(Hs(:,:,i), nrows, ncols);
    LHS(indi:indf, :) = real(Temp);
    RHS(indi:indf, :) = real(Hi) - Sys.D;
    LHS((step + indi):(step + indf), :) = imag(Temp);
    RHS((step + indi):(step + indf), :) = imag(Hi);
  end
  Sys.B = LHS \ RHS;
end
