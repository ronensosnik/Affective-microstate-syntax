function tf = lt_alpha(p, alpha)
% true if p is finite and < alpha
tf = ~isnan(p) && isfinite(p) && p < alpha;
end
