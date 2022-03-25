function [D, reHs] = desc_guessD(reHs, opts)

  % Usage : [D, reHs] = desc_guessD(real(Hs), opts)
  %
  % Guess a real D matrix from the given data.
  % The matrix is substracted from Hs and given
  % in reHs.
  %
  % opts.type is optional. It could be inf.
  
  if nargin == 1
    opts.type = 'norm';
    opts.isSymmetric = false;
  end
  if ~isfield(opts, 'type')
    opts.type = 'norm';
  end
  if ~isfield(opts, 'isSymmetric')
    opts.isSymmetric = false;
  end
  [nrows, ncols, nfreqs] = size(reHs);
  vs = zeros(nfreqs);
  if opts.type == 'norm'
    for idx = 1:nfreqs
      vs(idx) = norm(reHs(:,:,idx));
    end
    [~, idx] = sort(vs);
    D = reshape(reHs(:,:,idx(1)), nrows, ncols) + 0.01 * eye(nrows, ncols);
  elseif opts.type == 'cond'
    for idx = 1:nfreqs
      vs(idx) = cond(reHs(:,:,idx));
    end
    [~, idx] = sort(vs);
    D = reshape(reHs(:,:,idx(1)), nrows, ncols);
  end
  if opts.isSymmetric
    D = 0.5 * (D + transpose(D));
  end
  if norm(D) >= 1
    [U, S, V] = svd(D);
    sigmas = diag(S);
    allowed = ones(size(sigmas)) - eps;
    sigmas = min(allowed, sigmas);
    D = U * diag(sigmas) * V';
  end
  for idx = 1:nfreqs
    reHs(:,:,idx) = reHs(:,:,idx) - D;
  end
end
