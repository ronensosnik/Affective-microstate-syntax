function T = normalize_table(T)
% Ensure it's a table; convert categoricals to strings; logical to double.
if istable(T)
    % ok
elseif isnumeric(T)
    T = array2table(T);
else
    error('Expected table or numeric array.');
end

for j = 1:width(T)
    v = T.(j);
    if iscategorical(v)
        T.(j) = string(v);
    elseif isstring(v)
        % leave
    elseif iscellstr(v)
        % leave
    elseif islogical(v)
        T.(j) = double(v);
    else
        % numeric or other: leave as-is
    end
end
end