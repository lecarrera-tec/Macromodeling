function [DF, LM] = desc_lm_heur(freqs, Hs, opts)
  % Usage : DF = desc_lm_heur(freqs, Hs, opts)
  %
  % Construct the model [A, B, C, D, E].
  %
  % freqs: Frequencies (1 x M)
  %
  % Hs: Transfer functions (NO x NI x M), where NO 
  % and NI are the number of the output and input 
  % system ports, and M the number of frequency 
  % samples.
  %
  % opts.ktimes : integer
  % If not present, the default is:
  %    ceil ( sqrt( 2 * nop * nip / nfreqs ) )
  %    where nop, nip, nfreqs are the number of output 
  %    and input poles, and the number of frequencies
  %    respectively.
  % Size of the blocks. If -1 or inf, the maximum possible.
  % If not, should be a positive integer less or
  % equal than min(nip, nop).
  %
  % opts.npoles : integer
  % Number of poles (size of the system)
  % If not present or <= 0, uses the value of 
  % alpha to estimate it.
  %
  % opts.alpha : positive float (default : 1)
  % It works only if npoles <= 0.
  % alpha is the ratio of the total slope. If 1,
  % the same. If alpha < 1, bigger model.
  % If alpha > 1, smaller model (be careful).
  %
  % opts.vfit : boolean (default : false)
  % If true, then uses an even number of poles
  % for each column, so that VF could be made
  % exactly of the same size.

  global fid;

  if ~isfield(opts, 'dtype')
    opts.dtype = 'norm';
  end

  if ~isfield(opts, 'vfit')
    opts.vfit = false;
  end

  if ~isfield(opts, 'alpha')
    opts.alpha = 1;
  end

  DF.Notes = '';

  tic();
  reHs = real(Hs);
  imHs = imag(Hs);
  [nop, nip, nfreqs] = size(reHs);
  assert(nfreqs == size(freqs, 2));

  if ~isfield(opts, 'ktimes') || opts.ktimes < 1
    opts.ktimes = ceil(sqrt(2 * nop * nip / nfreqs));
  end
  opts.ktimes = min([opts.ktimes, nip, nop, floor(0.49 * nfreqs)]);
  ktimes = opts.ktimes;
  fprintf(fid, "ktimes = %d\n\n", ktimes);
    
  % 0. Extract D matrix.
  [DF.D, reHs] = desc_guessD(reHs, opts);

  kfreqs = ktimes * nfreqs;
  hkfreqs = 0.5 * kfreqs;

  R = spalloc(nip, kfreqs, hkfreqs);
  LM.W = zeros(nop, kfreqs);

  L = spalloc(kfreqs, nop, hkfreqs);
  LM.V = zeros(kfreqs, nip);

  II = eye(ktimes, ktimes);
  zz = zeros(ktimes, ktimes);
  omegas = 2 * pi * freqs;
  Lambda = kron(diag(omegas(1:2:end)), [zz, II; -II, zz]);
  M = kron(diag(omegas(2:2:end)), [zz, II; -II, zz]);

  % Constructing W, V, R and L matrices.
  Ii = eye(nip, nip);
  Io = eye(nop, nop);
  idx = (1:ktimes) - (ktimes + 1);
  iwv = (1:ktimes) - ktimes;
  for k = 1:2:nfreqs
    idx = idx + ktimes;
    iwv = iwv + ktimes;
    wdx = mod(idx,nip) + 1;
    vdx = mod(idx,nop) + 1;
    LM.W(:,iwv) = reHs(:,wdx,k);
    R(:,iwv) = Ii(:, wdx);
    LM.V(iwv,:) = reHs(vdx,:,k+1);
    L(iwv,:) = Io(vdx, :);
    iwv = iwv + ktimes;
    LM.W(:,iwv) =  imHs(:,wdx,k);
    LM.V(iwv,:) = -imHs(vdx,:,k+1);
  end; % for k = 1:2:nfreqs

  R = sqrt(2) * R;
  L = sqrt(2) * L;
  LM.W = sqrt(2) * LM.W;
  LM.V = sqrt(2) * LM.V;

  M2 = full(M * M);
  Lambda2 = full(Lambda * Lambda);
  T = LM.V * R - L * LM.W;

  % Solving the Sylvester equations "by hand".
  fprintf(fid, 'The sylvester equations ...\n');
  LM.IL = reshape(diag(M2), [], 1) - reshape(diag(Lambda2), 1, []);
  LM.sL = full(M*T*Lambda - L*LM.W*Lambda*Lambda + M2*LM.V*R) ./ LM.IL;
  LM.IL = full(M * T + T * Lambda) ./ LM.IL;

  T = omegas(1) * LM.IL - LM.sL;
  fprintf(fid, 'The SVD (%d x %d) ...\n', size(T,1), size(T,2));
  [LM.Y,LM.S,LM.X] = svd(T);

  sigmas = diag(LM.S);
  sigmas = sigmas / sigmas(1);
  if isfield(opts, 'npoles') && opts.npoles > 0
    n = opts.npoles;
  else
    n = sigor(sigmas, opts.alpha);
    plotesvd(sigmas, opts.prefix);
  end
  if opts.vfit
    n = nip * ceil(n / nip);
  end
  YY = LM.Y(:,1:n);
  XX = LM.X(:,1:n);
  DF.A = -YY' * LM.sL * XX; DF.B = YY' * LM.V;
  DF.C = LM.W * XX; DF.E = -YY' * LM.IL * XX;
  tt = toc();
  DF.Notes = [DF.Notes sprintf('Suggested system size = %d\n', n)];
  DF.Notes = [DF.Notes sprintf('DF Total Time = %g\n', tt)];
end % function desc_lm_heur


% alpha is the ratio of the total slope. If 1,
% the same. If alpha < 1, bigger model.
% If alpha > 1, smaller model (be careful).
function idx = sigor(sigmas, alpha)
  % plot(sigmas)
  global fid;
  sigmas = sigmas';
  n = length(sigmas);
  idx = 2:n;
  slope = 1 / (n - 1);
  deltas = sigmas(1:end-1) - sigmas(2:end);
  dslogic = deltas < alpha * slope;
  idx(dslogic) = [];
  outliers = 1.5 * iqr(idx) + quantile(idx, 0.75);
  fprintf(fid, "outliers = %.0f\n", outliers);
  idx = idx(idx < outliers);
  idx = idx(end);
end

function plotesvd(sigmas, fname)
global fid;
fprintf(fid, 'Graficando valores singulares ...\n');
fflush(fid);
data = [vec(1:length(sigmas)), vec(log10(sigmas))];
fname = ['~/NonReciprocal/SVDecomp/',  fname];
savename = [fname, '.data'];
save('-ascii', savename, 'data');
clf;
plot(log10(sigmas));
fprintf(fid, '... guardando grafica.\n\n');
fflush(fid);
savename = [fname, '.pdf'];
set(gcf, 'PaperSize', [10, 6])
saveas(gcf, savename);
end % function ploteigs
