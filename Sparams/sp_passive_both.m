function [DF, changed] = sp_passive_both(DF, omax=Inf, tol=1)
  % Usage: UDF = sp_passive_both(DF, omax=Inf, tol = 1)
  %
  % Updates DF.C matrix to construct a passive 
  % system. Hamiltonian perturbed algorithm.
  %
  % DF : DFtem constructed either with ss or dss.
  %
  % omax : maximum omega frequency of interest.
  % 
  % tol : Tolerance for the zero of the eigenvalues. 1
  % seems to be a good choice.
  %
  global fid;
  fprintf(fid, '@@ Full Global Passive for S-params @@\n');
  changed = false;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  opts.firsttime = true;
  opts.tol = 0;
  [omegas, sigmas] = sp_passv(DF, omax, opts)
  counter = 0;
  fprintf(fid, '** Hamiltonian search + local enforcement **\n');
  while ~isempty(omegas) && counter < 30
    changed = true;
    counter = counter + 1;
    fprintf(fid, 'counter = %d:', counter);
    DF = sp_lpe_list(DF, omegas, sigmas);
    [omegas, sigmas] = sp_passv(DF, omax);
    fprintf(fid, '\n');
  end
  if isempty(omegas)
    return;
  end
  fprintf(fid, '\n');

  %%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf(fid, '** Hamiltonian search + enforcement **\n');
  alpha = 0.3;
  
  Qcmt = DF.Qcmt;

  [nrows, ncols] = size(DF.D);

  SS = eye(nrows) - DF.D * DF.D';
  RR = eye(ncols) - DF.D' * DF.D;

  M11t = DF.B * DF.D';
  TT = SS \ DF.C;
  M11 = DF.A + M11t * TT;
  M21 = -DF.C' * TT;
  TT = RR \ DF.B';
  M12 = DF.B * TT;
  M22t = DF.D * TT;
  M22 = -DF.A' - DF.C' * M22t;
  MM = [M11, M12; M21, M22];
  if isdiag(DF.E) && prod(diag(DF.E)) == 1
    K0 = 1;
  else
    K0 = blkdiag(DF.E, DF.E');
  end


  fprintf(fid, ' @ Eigenvalues time ...\n');
  gg = tic();
  if K0 == 1
    [V, L] = eig(MM);
  else
    [V, L] = eig(MM, K0);
  end
  L = diag(L);
  fprintf(fid, "%f\n\n", toc(gg));

  % Forget the non-positive imag.
  indx = imag(L) > 0;
  L = L(indx);
  V = V(:,indx);

  % Only with almost zero real part.
  indx = find(abs(real(L)) < tol);
  L = imag(L(indx));
  V = V(:,indx);

  % We are not interested in higher freqs.
  indx = find(L <= omax);
  L = L(indx);
  V = V(:, indx);

  [L,indx] = sort(L);
  V = V(:,indx);
  n = length(L);
  fprintf(fid, '%d passivity violations ...\n', n);

  if n == 0
    fprintf(fid, 'The system is passive!!\n\n');
    return;
  end

  MP11t = -2 * DF.B * (SS \ DF.D);
  TT = RR \ DF.B';
  MP12 = (DF.B - 2 * TT') * TT;
  MP22t = 2 * (RR \ DF.D')' * TT;

  ssize = size(DF.A, 1);
  JJ = [zeros(ssize), -eye(ssize); eye(ssize), zeros(ssize)];

  while n > 0
    fprintf(fid, 'Freq. violations:');
    fprintf(fid, ' -- %.2f', L / (2e9 * pi));
    fprintf(fid, '\n');

    slopes = zeros(n,1);
    TT = SS \ DF.C;
    MP11 = MP11t * TT;
    MP21 = (-DF.C' - 2 * TT') * TT;
    MP22 = DF.C' * MP22t;
    MP = [MP11, MP12; MP21, MP22];
    for k = 1:n
      vk = V(:,k);
      v1 = vk(1:ssize);
      v2 = vk(ssize+1:end);
      slopes(k) = real(4 * imag(v1' * v2) / (vk' * JJ * MP * vk));
    end
    
    if slopes(n) > 0
      n = n - 1;
    end

    % We check first the first point.
    targets = zeros(n,1);
    if slopes(1) < 0
      targets(1) = 0;
    elseif n == 0 && omax < Inf
      targets = omax;
      n = 1;
      fprintf(fid, '\n @@ Una sola. Funcionara??? \n\n');
    elseif slopes(2) > 0
      targets(1) = L(2);
    else
      Delta = L(2) - L(1);
      targets(1) = L(1) + alpha * Delta;
    end

    % The rest of points. 
    for k = 2:n
      if slopes(k) > 0
        Delta = L(k+1) - L(k);
        temp = slopes(k+1);
      else
        Delta = L(k-1) - L(k);
        temp = slopes(k-1);
      end

      % If they have the same slope, jump to the
      % next point. If they have different slope,
      % look a point in the middle.
      if temp * slopes(k) > 0
        temp = 1;
      else
        temp = alpha;
      end
      targets(k) = L(k) + temp * Delta;
    end % for k = 2:n

    Zt = zeros(prod(size(DF.C)),n);
    eta = zeros(n,1);
    for k = 1:n
      v1 = V(1:ssize,k);
      v2 = V(ssize+1:end,k);
      yk = SS * DF.C * v1 + DF.D * RR * DF.B' * v2;
      Zt(:,k) = real(kron(Qcmt * v1, conj(yk)));
      eta(k) = (L(k) - targets(k)) * imag((DF.E * v1)' * v2);
    end

    [Q, R] = qr(Zt, 0);
    Chi = reshape(Q * (R' \ eta), size(DF.C));
    DF.C = DF.C + Chi * Qcmt;

    TT = SS \ DF.C;
    M11 = DF.A + M11t * TT;
    M21 = -DF.C' * TT;
    M22 = -DF.A' - DF.C' * M22t;
    MM = [M11, M12; M21, M22];

    if K0 == 1
      [V, L] = eig(MM);
    else
      [V, L] = eig(MM, K0);
    end
    L = diag(L);

    % Forget the non-positive imag.
    indx = imag(L) > 0;
    L = L(indx);
    V = V(:,indx);

    % Only with almost zero real part.
    indx = abs(real(L)) < tol;
    L = imag(L(indx));
    V = V(:,indx);

    % We are not interested in higher freqs.
    indx = find(L <= omax);
    L = L(indx);
    V = V(:, indx);

    [L,indx] = sort(L);
    V = V(:,indx);
    n = length(L);
    fprintf(fid, '%d passivity violations ...\n', n);
  end
  fprintf(fid, '\n');
end
