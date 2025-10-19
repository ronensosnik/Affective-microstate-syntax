function val = getfield_safe(S, f, def)
try
    if isfield(S,f) && ~isempty(S.(f)), val = double(S.(f)); else, val = def; end
catch, val = def;
end
end