function OUT = load_out_struct(matFile)

% Load the .mat and try to find a struct with fields detail & summary.

S = load(matFile);
OUT = [];
fn = fieldnames(S);

% Case 1: variables named 'detail' and 'summary'

if ismember('detail', fn) && istable(S.detail) && ismember('summary', fn) && istable(S.summary)
    OUT.detail = S.detail;
    OUT.summary = S.summary;
    return;
end

% Case 2: a single struct variable that contains fields detail/summary

for k = 1: numel(fn)
    val = S.(fn{k});

    if isstruct(val) && isfield(val, 'detail') && isfield(val, 'summary') && istable(val.detail) && istable(val.summary)
        OUT = val;
        return;
    end
end

% Case 3: a single variable that itself is the OUT struct (common name)
% Try known names like Coverage_Baseline_vs_meds, etc.

for k = 1: numel(fn)
    val = S.(fn{k});

    if isstruct(val) && isfieldNested(val, {'detail'}) && isfieldNested(val, {'summary'})
        try
            if istable(val.detail) && istable(val.summary)
                OUT = val;
                return;
            end
        catch
            % If stored as old-style, attempt to convert (rare)
        end
    end
end
end