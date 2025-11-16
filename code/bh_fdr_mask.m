function [sigMask, qMat] = bh_fdr_mask(pMat, alpha)

if nargin < 2
    alpha = 0.05;
end

v = pMat(:);
keep = ~isnan(v);
q = NaN(size(v));

if any(keep)
    q(keep) = bh_adjust(v(keep));
end

qMat = reshape(q, size(pMat));
sigMask = qMat < alpha;
end
