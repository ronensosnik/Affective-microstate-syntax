function v = safeget(S, varargin)
% Nested safe field access; returns NaN if missing; for tables returns [].
try
    for k = 1:numel(varargin)
        S = S.(varargin{k});
    end
    v = S;
catch
    if isstruct(S) || iscell(S)
        v = [];
    else
        v = NaN;
    end
end
end