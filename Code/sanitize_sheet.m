function s = sanitize_sheet(s)

% Excel sheet name must be <= 31 chars, no []:*?/\

s = regexprep(s, '[][\\/*?:]', '_');

if strlength(s) > 31
    s = extractBefore(s, 32);
end

end
