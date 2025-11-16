function C = dumpPairTable(T, measure, ERP, Micro, scope, factor, ctxVar, pipeline_label)

% Returns a cell array of rows for the output table.

C = {};

if isempty(T) || ~istable(T)
    return;
end

if ~ismember('pFDR', T.Properties.VariableNames)
    % If your saved tables lack pFDR (unlikely with your pipeline), skip.
    return;
end

sig = ~isnan(T.pFDR) & T.pFDR < 0.05;
T = T(sig, :);

if isempty(T)
    return;
end

% Determine effect columns by measure

switch lower(measure)
    case 'duration'
        eff = T.Delta;
        lo = T.Delta_Low;
        hi = T.Delta_High;
        units = "ms";

    case 'coverage'
        eff = T.OR;
        lo = T.OR_Low;
        hi = T.OR_High;
        units = "OR";

    otherwise  % 'occurrence'
        eff = T.RR;
        lo = T.RR_Low;
        hi = T.RR_High;
        units = "RR";
end

% Figure out level columns and optional context

if ismember('Level1', T.Properties.VariableNames)

    L1 = string(T.Level1);
    L2 = string(T.Level2);
    ctx = repmat("", height(T), 1);

elseif strcmpi(factor, 'Group') && all(ismember({'Group_1', 'Group_2'}, T.Properties.VariableNames))

    L1 = string(T.Group_1);
    L2 = string(T.Group_2);

    if ~isempty(ctxVar) && ismember(ctxVar, T.Properties.VariableNames)
        ctx = string(T.(ctxVar));
    else
        ctx = repmat("", height(T), 1);
    end

elseif strcmpi(factor, 'Condition') && all(ismember({'Condition_1', 'Condition_2'}, T.Properties.VariableNames))

    L1 = string(T.Condition_1);
    L2 = string(T.Condition_2);

    if ~isempty(ctxVar) && ismember(ctxVar, T.Properties.VariableNames)
        ctx = string(T.(ctxVar));
    else
        ctx = repmat("", height(T), 1);
    end
else

    % Fallback: attempt generic *_1 / *_2 detection

    cols = T.Properties.VariableNames;
    lefts = cols(endsWith(cols, '_1'));
    rights = cols(endsWith(cols, '_2'));

    if ~isempty(lefts) && ~isempty(rights)
        L1 = string(T.(lefts{1}));
        L2 = string(T.(rights{1}));
    else
        L1 = repmat("", height(T), 1);
        L2 = L1;
    end
    ctx = repmat("", height(T), 1);
end

% Build rows

C = cell(height(T), 15);

for i = 1: height(T)
    C(i, :) = {char(measure), char(ERP), char(Micro), char(scope), char(factor), char(ctx(i)), char(L1(i)), char(L2(i)), eff(i), lo(i), hi(i), T.pValue(i), T.pFDR(i), char(units), char(pipeline_label)};
end
end