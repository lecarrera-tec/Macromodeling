function [A, bs, x0] = feasible(A, bs)
% Usage: [nA, nbs, x0] = feasible(A, bs)
% Find x0 such that nA * x0 >= nbs
% nA and nbs are reordered versions of A y bs
% in terms of bs values.
  n = length(bs);
  [xs, ixs] = sort(bs);
  ixs = ixs(n:-1:1);
  bs = bs(ixs);
  A  = A(ixs,:);
  value = bs(1) / sum(abs(A(1,:)))
  x0 = value * sign(A(1,:))';
  assert(A(1,:) * x0 >= bs(1));

  i = 2;
  eqsigns = sign(A(1,:)) == sign(A(2,:));
  feasible = sum(eqsigns);
  while (feasible > 0 && i < n)
    missed = bs(i) - A(i,:) * x0;
    if (missed > 0)
      value = missed / sum(abs(eqsigns .* A(i,:)))
      x0 = x0 + (eqsigns .* (value .* sign(A(i,:))))';
      %A(1:i,:) * x0 - bs(1:i,:)
      %assert(A(1:i,:) * x0 >= bs(1:i,:));
    end
    i = i + 1;
    eqsigns = eqsigns .* (sign(A(i,:)) == sign(A(1,:)));
    feasible = sum(eqsigns);
  end
  A = A(1:i-1,:);
  bs = bs(1:i-1,:);
end
