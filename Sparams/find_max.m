function [xmax, ymax] = find_max(xs, ys, Sys)
  % Usage: [xmax, ymax] = find_max(xs, ys, Sys)
  % Finds the frequency for the maximum norm of the
  % linear system Sys.
  %
  % xs : 3 frequencies
  % ys : already computed norms for the given frequencies
  % Sys : linear system.
  %
  % As the difference for the values is minimum, we
  % will use SVD to solve the system with greater
  % precision.

  % We will solve a quadratic equation.
  M = zeros(3,3);
  tol = 1e-3;
  delta = inf;

  % Choose the maximum norm between
  % the given values.
  [ymax, ind] = max(ys);
  xmax = xs(ind);

  while (delta > tol)
    % Constructing the quadratic equation.
    % cs are the coefficients
    M(:,1) = reshape(xs .* xs, [], 1);
    M(:,2) = reshape(xs, [], 1);
    M(:,3) = ones(3,1);
    [U,S,V] = svd(M);
    cs = V * inv(S) * U' * reshape(ys, [], 1);

    % Is a convex (concave upwards) function.
    % Take the maximum and leave.
    if cs(1) >= 0
      [yi,ind] = max(ys);
      if ymax < yi
        xmax = xs(ind);
        ymax = yi;
      end
      return;
    end

    % the x coordinate of the vertex.
    h = -0.5 * cs(2) / cs(1);

    % If the omega is outside the
    % given range, find maximum and leave.
    if h < xs(1) || h > xs(3)
      [yi,ind] = max(ys);
      if ymax < yi
        xmax = xs(ind);
        ymax = yi;
      end
      return;
    end

    % Getting the norm.
    Hi = evalfr(Sys, 1j * h);
    k = norm(Hi);
    if ymax < k
      xmax = h;
      ymax = k;
    end

    % Distance in x from the central point to the vertex.
    delta = abs(xs(2) - h) / h;

    % Constructing a smaller interval.
    if h < xs(2)
      xs(3) = xs(2);
      ys(3) = ys(2);
    else
      xs(1) = xs(2);
      ys(1) = ys(2);
    end
    xs(2) = h;
    ys(2) = k;
  end
end
