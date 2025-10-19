function S = nested_get(S, keys)
for kk = 1:numel(keys)
    k = keys{kk};
    if ~isfield(S,k), S = []; return; end
    S = S.(k); if isempty(S), return; end
end
end