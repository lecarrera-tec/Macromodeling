% 1. Find x suboptimal such that \sum(A * x) = \sum(bs)
% bs > 0;
% 2. Find the common signs, and then
% start to fill xs from the bigger bs.
function xs = mixed(A, bs)
  [nrows, ncols] = size(A);

  % Trivial case.
  if nrows == 1
    denom = A * A';
    xs = A' * (bs(1) / denom);
    return
  end

  % First part. Quick sub-filling.
  siga = sum(A);
  sigb = sum(bs);
  deno = siga * siga';
  xs = siga' * (sigb / deno);

  % Update the necesities.
  bs = bs - A * xs;

  % Second part. Complete filling.
  % Find common signs.
  i = 2;
  oldeqsg = ones(1, ncols);
  eqsigns = sign(A(1,:)) == sign(A(2,:));
  feasible = sum(eqsigns);

  while (2 * feasible > ncols && i < nrows)
    i = i + 1;
    oldeqsg = eqsigns;
    eqsigns = oldeqsg .* (sign(A(i,:)) == sign(A(1,:)));
    feasible = sum(eqsigns);
  end

  if (2 * feasible < ncols)
    eqsigns = oldeqsg;
    nrows = i - 1;
    A = A(1:nrows,:);
    bs = bs(1:nrows,:);
  end
  clear oldeqsg nrows

  % Initialization
  [bi,ii] = max(bs);

  % While it has sense.
  while (size(bs,1) > 0)
    ai = eqsigns .* A(ii,:);
    deno = ai * ai';
    xupd = ai * ((bi * (1.0 + 2.0 * eps)) / deno);
    xs = xs + xupd';

    % Update the necesities.
    bs = bs - A * xupd';

    % Remove the unnecessary rows.
    A(bs<=eps,:) = [];
    bs(bs<=eps,:) = [];

    % Update the info.
    [bi,ii] = max(bs);
  end

end
