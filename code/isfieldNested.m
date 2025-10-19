function tf = isfieldNested(S, pathCell)
tf = true;
try
    for i = 1:numel(pathCell)
        S = S.(pathCell{i});
    end
catch
    tf = false;
end
end
