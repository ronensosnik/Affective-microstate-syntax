function s = pct(a,b)

if b <= 0
    s = "NA";
else
    s = sprintf('%.1f%%', 100*a/b);
end
end
