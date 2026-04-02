function lines = structTreeSummary(S, name)
% structTreeSummary : Recursively summarize struct fields (tree view).
%
% Usage:
%   structTreeSummary(S)
%   structTreeSummary(S, 'S')
%   lines = structTreeSummary(S, 'S')   % returns cellstr lines

if nargin < 2
    name = inputname(1);
    if isempty(name)
        name = 'S';
    end
end

lines = {};
walk(S, name, 0);

if nargout == 0
    fprintf('%s\n', lines{:});
end

    function walk(x, path, indent)
        prefix = [repmat('  ', 1, indent) '- '];
        lines{end + 1} = [prefix path ' : ' describeValue(x)];

        if isstruct(x)
            if numel(x) ~= 1
                % For non-scalar struct arrays, don't descend into every element by default.
                % Show fields once, and optionally descend into the first element.
                f = fieldnames(x);
                for i = 1 : numel(f)
                    v1 = x(1).(f{i});
                    walk(v1, [path '(' num2str(size(x, 1)) 'x' num2str(size(x, 2)) ').' f{i}], indent + 1);
                end
                return
            end

            f = fieldnames(x);
            for i = 1 : numel(f)
                v = x.(f{i});
                walk(v, [path '.' f{i}], indent + 1);
            end

        elseif iscell(x)
            % Optionally show a small preview of cell contents types (first few).
            nPreview = min(5, numel(x));
            for k = 1 : nPreview
                walk(x{k}, [path '{' num2str(k) '}'], indent + 1);
            end
            if numel(x) > nPreview
                lines{end + 1} = [repmat('  ', 1, indent + 1) '- ' path '{...} : (cell preview truncated)'];
            end
        end
    end

    function s = describeValue(x)
        c = class(x);
        sz = size(x);
        szStr = [num2str(sz(1)) 'x' num2str(sz(2))];
        if numel(sz) > 2
            for j = 3 : numel(sz)
                szStr = [szStr 'x' num2str(sz(j))];
            end
        end

        if isstruct(x)
            s = ['struct ' szStr];
        elseif iscell(x)
            s = ['cell ' szStr];
        elseif ischar(x)
            s = ['char ' szStr];
        elseif isstring(x)
            s = ['string ' szStr];
        elseif isa(x, 'table')
            s = ['table ' szStr];
        elseif isa(x, 'timetable')
            s = ['timetable ' szStr];
        else
            s = [c ' ' szStr];
        end
    end
end