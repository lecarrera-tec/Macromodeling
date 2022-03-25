function Sys = desc_updateC(Sys, Hs, omegas)
  % Usage: Sys = desc_updateC(Sys, Hs, omegas)
  % Update C to minimize the error.
  % It is an overdetermined system of the form
  % C * LHS = RHS \equiv LHS^T * C^T = RHS^T

  % Construct first matrices LHS and RHS.
  global fid;
  fprintf(fid, 'Updating C ...\n');

  nfreqs = size(omegas, 2);
  ssize = size(Sys.A, 1);
  [nrows, ncols] = size(Sys.D);
  step = ncols * nfreqs;
  SysPsy = Sys;
  SysPsy.C = eye(ssize);
  SysPsy.D = 0;
  Psy = desc_eval(omegas, SysPsy);

  indi = 0;
  indf = 0;
  LHS = zeros(ssize, 2 * step);
  RHS = zeros(nrows, 2 * step);
  for i = 1:nfreqs
    indi = indf + 1;
    indf = indf + ncols;
    Temp = reshape(Psy(:, :, i), ssize, ncols);
    Hi = reshape(Hs(:,:,i), nrows, ncols) - Sys.D;
    LHS(:, indi:indf) = real(Temp);
    RHS(:, indi:indf) = real(Hi);
    LHS(:, (step + indi):(step + indf)) = imag(Temp);
    RHS(:, (step + indi):(step + indf)) = imag(Hi);
  end
  Sys.C = RHS / LHS;
end
