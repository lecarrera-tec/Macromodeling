function DF = sp_lpe_list(DF, omegas, sigmas, delta)

% Usage: Usys = sp_lpe_list(DF, omegas, sigmas, delta = 5e-12)
%
% Computes the delta DF.C matrix. Local passivity algorithm for the given omegas.
%
% DF : Descriptor Form or Space State system, with fields 
% DF.A, DF.B, DF.C, DF.D, (possibly DF.E) and Qcmt.
% omegas : Frequencies with passivity violations.
% sigmas : Norms for each given frequency.
% delta  : Value below 1. Default is 5e-12.

global fid;

if nargin < 4
  delta = 5e-12;
end

ssize = size(DF.A, 1);
[nports, ncols] = size(DF.D);
assert(nports == ncols);
PhiDF = DF;
PhiDF.C = eye(ssize);
PhiDF.D = 0;
Qcmt = DF.Qcmt;
hinf = inf;

upward = 0;
while ~isempty(omegas) && upward < 10
  last_hinf = hinf;
  hinf = max(sigmas) - 1;
  fprintf(fid, ' %.5e', hinf);
  if (hinf >= last_hinf)
    upward = upward + 1;
    fprintf(fid, ' (upward = %d)', upward);
  end
  last_hinf = hinf;
  lo = length(omegas);
  if lo > 19
    idx = [sigmas(1:end-1) >= sigmas(2:end), 1];
    idx = idx .* [1, sigmas(2:end) >= sigmas(1:end-1)];
    if nnz(idx) < 6
      idx = idx + [0, idx(1:end-1)] + [idx(2:end), 0];
    end
    if nnz(idx) < 10
      temp = 1 + find(idx(2:end-1));
      for ii = temp
        [~, jj] = max([sigmas(ii - 1), 0, sigmas(ii + 1)]);
        idx(ii + jj - 2) = 1;
      end
      if idx(1) > 0
        idx(2) = 1;
      end
      if idx(end) > 0
        idx(end-1) = 1;
      end
    end
    idx = logical(idx);
    subset = omegas(idx);
    sigmas = sigmas(idx);
  else
    subset = omegas;
  end
  
  n = length(subset);
  Phi = desc_eval(subset, PhiDF);
  Z = zeros(n, ssize * nports);
  bb = zeros(n, 1);
  for kk = 1:n
    Phii = reshape(Phi(:, :, kk), ssize, nports);
    Hi = DF.C * Phii + DF.D;
    [U, S, V] = svd(Hi);
    bb(kk) = 1 - delta - S(1);
    Left = Qcmt * Phi(:,:,kk);
    xji = real(kron(Left * V(:,1), conj(U(:,1))));
    Z(kk,:) = xji';
  end

  Xi = mixed(-Z, -bb);
  DF.C = DF.C + reshape(Xi, size(DF.C)) * Qcmt;

  sigmas = desc_norms(DF, omegas);
  idx = sigmas >= 1;
  sigmas = sigmas(idx);
  omegas = omegas(idx);
end

end
