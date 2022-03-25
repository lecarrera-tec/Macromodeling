function [xs, ys] = del_repeated(xs, ys)
% Usage [nxs, nys] = del_repeated(xs, ys)
% xs is an ordered list. ys is related to xs.
% It search for repeated values in xs, and
% removes the respective entries for both lists.

zs = [];
ns = length(xs);
for ii = 2:ns
  if xs(ii-1) == xs(ii)
    zs(end+1) = ii;
  end
end
xs(zs) = [];
ys(zs) = [];

end
