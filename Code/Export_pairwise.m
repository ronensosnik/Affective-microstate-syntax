function Export_pairwise(OUT, measure, pipeline_label, out_csv)

% EXPORT_PAIRWISE
% Export ONLY pairwise comparisons that passed BOTH:
%   (1) ERP-wise BH-FDR omnibus significance, and
%   (2) within-family (cell) BH-FDR significance for the pairwise contrast.
%
% Output columns:
% Measure, ERP, Microstate, Scope, Factor, Context, Level1, Level2, Effect,
% CI_low, CI_high, p, pFDR, Units, Pipeline
%
% Formatting:
%   Effect / CI columns: 3 decimals (e.g., 0.524)
%   p and pFDR: <0.001 if smaller; else 3 decimals

if ischar(OUT) || isstring(OUT)
    OUT = load_out_struct_from_mat(OUT);
end

measure = validatestring(measure, {'Duration', 'Coverage', 'Occurrence', 'Num_occurrence'});

[measureField, measureLabel, unitsLabel] = resolve_measure_field_and_units(OUT, measure);

if ~isfield(OUT, 'info') || ~isfield(OUT.info, 'ERPs') || ~isfield(OUT.info, 'Microstates')
    error('OUT must contain OUT.info.ERPs and OUT.info.Microstates.');
end

if isfield(OUT.info, 'alphaFDR') && isfinite(OUT.info.alphaFDR)
    alphaFDR = OUT.info.alphaFDR;
else
    alphaFDR = 0.05;
end

ERPs = OUT.info.ERPs;
Micros = OUT.info.Microstates;

if ~isfield(OUT, measureField)
    error('OUT does not contain field "%s".', measureField);
end

rows = struct('Measure', {}, 'ERP', {}, 'Microstate', {}, 'Scope', {}, 'Factor', {}, ...
              'Context', {}, 'Level1', {}, 'Level2', {}, 'Effect', {}, 'CI_low', {}, ...
              'CI_high', {}, 'p', {}, 'pFDR', {}, 'Units', {}, 'Pipeline', {});

for e = 1: numel(ERPs)
    ERPname = ERPs{e};

    for m = 1: numel(Micros)
        Microname = Micros{m};

        R = OUT.(measureField){e, m};

        if isempty(R) || ~isstruct(R) || ~isfield(R, 'pairwise') || isempty(R.pairwise)
            continue;
        end

        pInt = get_nested(R, {'tests', 'interaction', 'pFDR'});
        pGrp = get_nested(R, {'tests', 'group', 'pFDR'});
        pCond = get_nested(R, {'tests', 'condition', 'pFDR'});

        if isfinite(pInt) && pInt < alphaFDR

            if isfield(R.pairwise, 'Group_within_Condition')
                P = R.pairwise.Group_within_Condition;
                rows = append_pairwise_rows_simple_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, ...
                                                      'Group within Condition', 'Group', 'Condition');
            end

            if isfield(R.pairwise, 'Condition_within_Group')
                P = R.pairwise.Condition_within_Group;
                rows = append_pairwise_rows_simple_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, ...
                                                      'Condition within Group', 'Condition', 'Group');
            end

        else

            if isfinite(pGrp) && pGrp < alphaFDR && isfield(R.pairwise, 'Group')
                P = R.pairwise.Group;
                rows = append_pairwise_rows_main_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, ...
                                                     'Group (marginal)', 'Group', string(missing));
            end

            if isfinite(pCond) && pCond < alphaFDR && isfield(R.pairwise, 'Condition')
                P = R.pairwise.Condition;
                rows = append_pairwise_rows_main_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, ...
                                                     'Condition (marginal)', 'Condition', string(missing));
            end

        end

    end
end

if isempty(rows)
    warning('No pairwise results met BOTH ERP-FDR and within-family FDR for measure "%s".', measureLabel);
    T = table();
else
    T = struct2table(rows);

    T.Measure = string(T.Measure);
    T.ERP = string(T.ERP);
    T.Microstate = string(T.Microstate);
    T.Scope = string(T.Scope);
    T.Factor = string(T.Factor);
    T.Context = string(T.Context);
    T.Level1 = string(T.Level1);
    T.Level2 = string(T.Level2);
    T.Units = string(T.Units);
    T.Pipeline = string(T.Pipeline);

    % Sort while numeric values are still numeric
    T = sortrows(T, {'Measure', 'ERP', 'Microstate', 'Scope', 'Factor', 'Context', 'pFDR', 'p'});

    % Format numeric columns to strings (3 decimals; p/pFDR special rule)
    T.Effect = fmt3(T.Effect);
    T.CI_low = fmt3(T.CI_low);
    T.CI_high = fmt3(T.CI_high);
    T.p = fmt_p(T.p);
    T.pFDR = fmt_p(T.pFDR);
end

[~, ~, ext] = fileparts(out_csv);

if any(strcmpi(ext, {'.xlsx', '.xls', '.xlsm'}))
    writetable(T, out_csv, 'Sheet', 'Pairwise_sig');
else
    writetable(T, out_csv);
end

end


function OUTs = load_out_struct_from_mat(matFile)

S = load(matFile);
fn = fieldnames(S);

if isempty(fn)
    error('MAT file "%s" contains no variables.', matFile);
end

OUTs = [];

for k = 1: numel(fn)
    v = S.(fn{k});
    if isstruct(v) && isfield(v, 'info')
        OUTs = v;
        return;
    end
end

v = S.(fn{1});
if isstruct(v)
    OUTs = v;
    return;
end

error('Could not find an OUT struct in "%s".', matFile);

end


function [measureField, measureLabel, unitsLabel] = resolve_measure_field_and_units(OUT, measureIn)

measureField = measureIn;
measureLabel = measureIn;

if strcmpi(measureIn, 'Occurrence')
    if isfield(OUT, 'Num_occurrence')
        measureField = 'Num_occurrence';
    else
        measureField = 'Occurrence';
    end
    measureLabel = 'Occurrence';
end

if strcmpi(measureIn, 'Num_occurrence')
    measureField = 'Num_occurrence';
    measureLabel = 'Occurrence';
end

if strcmpi(measureIn, 'Coverage')
    unitsLabel = 'OR';
elseif strcmpi(measureLabel, 'Occurrence')
    unitsLabel = 'RR';
else
    unitsLabel = 'ms';
end

end


function val = get_nested(S, path)

val = NaN;

try
    for i = 1: numel(path)
        S = S.(path{i});
    end
    val = S;
catch
end

end


function rows = append_pairwise_rows_main_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, scopeStr, factorStr, contextVal)

if isempty(P) || ~istable(P)
    return;
end

[pRaw, pFDR] = get_p_columns(P);

mask = isfinite(pFDR) & (pFDR < alphaFDR);
if ~any(mask)
    return;
end

est = get_col_or_nan(P, 'Estimate');
se = get_col_or_nan(P, 'SE');

idxKeep = find(mask);

for k = 1: numel(idxKeep)

    i = idxKeep(k);

    contrastStr = get_contrast_string(P, i, factorStr);
    [L1, L2] = split_contrast_safe(contrastStr);

    [eff, lo, hi] = effect_and_ci(est(i), se(i), measureLabel);

    r = struct();
    r.Measure = measureLabel;
    r.ERP = ERPname;
    r.Microstate = Microname;
    r.Scope = scopeStr;
    r.Factor = factorStr;
    r.Context = contextVal;
    r.Level1 = L1;
    r.Level2 = L2;
    r.Effect = eff;
    r.CI_low = lo;
    r.CI_high = hi;
    r.p = pRaw(i);
    r.pFDR = pFDR(i);
    r.Units = unitsLabel;
    r.Pipeline = pipeline_label;

    rows(end + 1, 1) = r;
end

end


function rows = append_pairwise_rows_simple_sig(rows, P, alphaFDR, measureLabel, unitsLabel, pipeline_label, ERPname, Microname, scopeStr, withinFactor, byFactor)

if isempty(P) || ~istable(P)
    return;
end

[pRaw, pFDR] = get_p_columns(P);

mask = isfinite(pFDR) & (pFDR < alphaFDR);
if ~any(mask)
    return;
end

est = get_col_or_nan(P, 'Estimate');
se = get_col_or_nan(P, 'SE');

if ~any(strcmp(P.Properties.VariableNames, byFactor))
    error('Pairwise table missing context factor "%s".', byFactor);
end

idxKeep = find(mask);

for k = 1: numel(idxKeep)

    i = idxKeep(k);

    ctx = string(P.(byFactor)(i));
    contrastStr = get_contrast_string(P, i, withinFactor);
    [L1, L2] = split_contrast_safe(contrastStr);

    [eff, lo, hi] = effect_and_ci(est(i), se(i), measureLabel);

    r = struct();
    r.Measure = measureLabel;
    r.ERP = ERPname;
    r.Microstate = Microname;
    r.Scope = scopeStr;
    r.Factor = withinFactor;
    r.Context = ctx;
    r.Level1 = L1;
    r.Level2 = L2;
    r.Effect = eff;
    r.CI_low = lo;
    r.CI_high = hi;
    r.p = pRaw(i);
    r.pFDR = pFDR(i);
    r.Units = unitsLabel;
    r.Pipeline = pipeline_label;

    rows(end + 1, 1) = r;
end

end


function x = get_col_or_nan(T, colName)

if any(strcmp(T.Properties.VariableNames, colName))
    x = T.(colName);
else
    x = nan(height(T), 1);
end

end


function [pRaw, pFDR] = get_p_columns(P)

pRaw = nan(height(P), 1);
pFDR = nan(height(P), 1);

if any(strcmp(P.Properties.VariableNames, 'pValue'))
    pRaw = P.pValue;
elseif any(strcmp(P.Properties.VariableNames, 'p'))
    pRaw = P.p;
end

if any(strcmp(P.Properties.VariableNames, 'pFDR'))
    pFDR = P.pFDR;
elseif any(strcmp(P.Properties.VariableNames, 'qValue'))
    pFDR = P.qValue;
elseif any(strcmp(P.Properties.VariableNames, 'pAdj'))
    pFDR = P.pAdj;
end

end


function s = get_contrast_string(P, i, factorStr)

s = "";

if any(strcmp(P.Properties.VariableNames, factorStr))
    s = string(P.(factorStr)(i));
    return;
end

if any(strcmp(P.Properties.VariableNames, 'Contrast'))
    s = string(P.Contrast(i));
    return;
end

for j = 1: width(P)
    if isstring(P{:, j}) || iscategorical(P{:, j})
        s = string(P{i, j});
        return;
    end
end

end


function [L1, L2] = split_contrast_safe(s)

s = strtrim(string(s));

s = replace(s, char(8211), '-');
s = replace(s, char(8212), '-');

parts = split(s, " - ");

if numel(parts) >= 2
    L1 = strtrim(parts(1));
    L2 = strtrim(parts(2));
else
    L1 = s;
    L2 = "";
end

end


function [eff, lo, hi] = effect_and_ci(est, se, measureLabel)

if ~isfinite(est) || ~isfinite(se)
    eff = NaN;
    lo = NaN;
    hi = NaN;
    return;
end

crit = 1.96;

if strcmpi(measureLabel, 'Duration')
    eff = est;
    lo = est - crit * se;
    hi = est + crit * se;
else
    eff = exp(est);
    lo = exp(est - crit * se);
    hi = exp(est + crit * se);
end

end


function s = fmt3(x)

% Format numeric vector as strings with 3 decimals; NaN -> "NA"
s = strings(size(x));

for i = 1: numel(x)
    if isfinite(x(i))
        s(i) = sprintf('%.3f', x(i));
    else
        s(i) = "NA";
    end
end

end


function s = fmt_p(x)

% Format p-values: <0.001 if smaller; else 3 decimals; NaN -> "NA"
s = strings(size(x));

for i = 1: numel(x)
    if ~isfinite(x(i))
        s(i) = "NA";
    elseif x(i) < 0.001
        s(i) = "<0.001";
    else
        s(i) = sprintf('%.3f', x(i));
    end
end

end
