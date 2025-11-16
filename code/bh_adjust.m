function q = bh_adjust(p, dim, method)

% bh_adjust  Benjamini–Hochberg (BH) / Benjamini–Yekutieli (BY) FDR q-values.
%   q = bh_adjust(p) % along first non-singleton dim, BH
%   q = bh_adjust(p, dim) % explicit dimension
%   q = bh_adjust(p, dim, method) % method: 'bh' (default) or 'by'
%
% INPUT
%   p: p-values in [0,1], any shape. NaNs are preserved.
%   dim: dimension to operate along (default: first non-singleton).
%   method: 'bh' (Benjamini–Hochberg, step-up) or 'by' (Benjamini–Yekutieli).
%
% OUTPUT
%   q: BH/BY adjusted q-values, same size as p, NaNs where p is NaN.
%
% NOTES
% - Monotonicity is enforced in the sorted domain, then mapped back.
% - Values are clipped to [0,1].
%
% REF
%  Benjamini & Hochberg (1995), Benjamini & Yekutieli (2001)

if nargin < 2 || isempty(dim)

    dim = find(size(p) ~= 1, 1);

    if isempty(dim)
        dim = 1;
    end
end

if nargin < 3 || isempty(method)
    method = 'bh';
end

method = validatestring(lower(method), {'bh',' by'});

% Work in double, keep original shape

sz = size(p);
q = NaN(sz, 'like', p);
x = double(p);

% Permute so that DIM is first, then collapse the rest into columns

order = [dim, setdiff(1:ndims(x), dim)];
xperm = permute(x, order);
n = size(xperm, 1);
m = prod(size(xperm)) / n;
X = reshape(xperm, n, m);

% Precompute BY constant if needed

if strcmp(method,'by')
    c_m = sum(1./(1:n));
else
    c_m = 1;
end

% Process each column independently

Q = NaN(size(X));

for j = 1: m
    pj = X(:, j);
    isn = isnan(pj);
    keep_idx = find(~isn);
    k = numel(keep_idx);

    if k == 0
        continue;
    end

    % Clamp to [0,1] just in case

    pj(keep_idx) = max(0, min(1, pj(keep_idx)));

    % Sort the kept p-values and remember their original positions

    [ps, sort_idx] = sort(pj(keep_idx), 'ascend');
    ranks = (1: k)';

    % Raw adjusted values

    adj = (k ./ ranks) .* ps * c_m;

    % Enforce monotonicity (non-increasing when traversing backwards)

    for ii = k - 1: -1: 1
        adj(ii) = min(adj(ii), adj(ii + 1));
    end

    % Clip to [0,1]
    adj = max(0, min(1, adj));

    % Map back to original row indices

    col_q = NaN(n, 1);
    orig_positions = keep_idx(sort_idx);
    col_q(orig_positions) = adj;

    Q(:, j) = col_q;
end

% Reshape back and unpermute

qperm = reshape(Q, size(xperm));
q = ipermute(qperm, order);

% Cast back like input (preserve single if provided)

if ~isa(q, class(p))
    q = cast(q, 'like', p);
end
end
