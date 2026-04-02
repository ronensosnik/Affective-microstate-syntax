function OUT = summarize_metric_sensitivity_robustness(Out_base, Out_sensitivity, measure)

% SUMMARIZE_METRIC_SENSITIVITY_ROBUSTNESS
% Compare baseline vs sensitivity LME results for a static metric.
%
% Adds robustness metrics:
%   - Omnibus agreement: percent agreement + Cohen's kappa (Interaction/Group/Condition)
%   - Pairwise robustness (when both pipelines have relevant pairwise tables):
%       * Jaccard overlap of significant pair sets (pFDR < alpha)
%       * Spearman correlation of effect sizes (Estimate; model scale)
%       * Median absolute effect-size change (|ΔEstimate|; model scale)
%   - EMM robustness (when both pipelines have relevant EMM tables):
%       * RMSE on model scale (Estimate)
%       * RMSE on back-transformed scale (Estimate_BT), when available
%
% Usage:
%   OUT = summarize_metric_sensitivity_robustness(Out_Duration_base, Out_Duration_sensitivity_smoothness, 'Duration');
%   OUT = summarize_metric_sensitivity_robustness(Out_Occurrence_base, Out_Occurrence_sensitivity_smoothness, 'Occurrence');
%   OUT = summarize_metric_sensitivity_robustness(Out_Coverage_base, Out_Coverage_sensitivity_smoothness, 'Coverage');

measure = validatestring(measure, {'Duration', 'Coverage', 'Occurrence', 'Num_occurrence'});

measureField = resolve_measure_field(Out_base, Out_sensitivity, measure);

ERPs = Out_sensitivity.info.ERPs;
Micros = Out_sensitivity.info.Microstates;

alpha = pick_alpha(Out_base, Out_sensitivity);

nCells = numel(ERPs) * numel(Micros);
rows(nCells, 1) = struct();

r = 0;

for e = 1: numel(ERPs)
    for m = 1: numel(Micros)

        r = r + 1;

        Rb = Out_base.(measureField){e, m};
        Rs = Out_sensitivity.(measureField){e, m};

        if isempty(Rb)
            Rb = struct();
        end

        if isempty(Rs)
            Rs = struct();
        end

        % ---------------- omnibus pFDRs ----------------

        pB_I = safeget(Rb, NaN, 'tests', 'interaction', 'pFDR');
        pB_G = safeget(Rb, NaN, 'tests', 'group', 'pFDR');
        pB_C = safeget(Rb, NaN, 'tests', 'condition', 'pFDR');

        pS_I = safeget(Rs, NaN, 'tests', 'interaction', 'pFDR');
        pS_G = safeget(Rs, NaN, 'tests', 'group', 'pFDR');
        pS_C = safeget(Rs, NaN, 'tests', 'condition', 'pFDR');

        bI = lt_alpha(pB_I, alpha);
        bG = lt_alpha(pB_G, alpha);
        bC = lt_alpha(pB_C, alpha);

        sI = lt_alpha(pS_I, alpha);
        sG = lt_alpha(pS_G, alpha);
        sC = lt_alpha(pS_C, alpha);

        lostI = bI & ~sI;
        gainedI = ~bI & sI;

        lostG = bG & ~sG;
        gainedG = ~bG & sG;

        lostC = bC & ~sC;
        gainedC = ~bC & sC;

        anyFDRchange = (lostI | gainedI) | (lostG | gainedG) | (lostC | gainedC);

        % ---------------- direction change (ordering) ----------------

        dirChanged = false;

        if bI && sI
            EB = safeget(Rb, [], 'emm', 'Group_by_Condition');
            ES = safeget(Rs, [], 'emm', 'Group_by_Condition');
            dirChanged = ~same_ordering_2F(EB, ES, "Group", "Condition");

        elseif ~bI && ~sI && bG && sG
            EB = safeget(Rb, [], 'emm', 'Group');
            ES = safeget(Rs, [], 'emm', 'Group');
            dirChanged = ~same_ordering_1F(EB, ES, "Group");

        elseif ~bI && ~sI && bC && sC
            EB = safeget(Rb, [], 'emm', 'Condition');
            ES = safeget(Rs, [], 'emm', 'Condition');
            dirChanged = ~same_ordering_1F(EB, ES, "Condition");
        end

        % ---------------- pairwise concordance metrics ----------------

        [jacc, rho, mad, famPW, nMatchPW] = pairwise_concordance(Rb, Rs, alpha);

        % ---------------- EMM RMSE metrics ----------------

        [rmseT, rmseBT, famEMM, nMatchEMM] = emm_rmse_concordance(Rb, Rs);

        % ---------------- store row ----------------

        rows(r).ERP = string(ERPs{e});
        rows(r).Map = string(Micros{m});
        rows(r).Measure = string(measureField);
        rows(r).AlphaFDR = alpha;

        rows(r).pFDR_Base_Int = pB_I;
        rows(r).pFDR_Base_Grp = pB_G;
        rows(r).pFDR_Base_Cond = pB_C;

        rows(r).pFDR_Sens_Int = pS_I;
        rows(r).pFDR_Sens_Grp = pS_G;
        rows(r).pFDR_Sens_Cond = pS_C;

        rows(r).Base_Int = bI;
        rows(r).Sensitivity_Int = sI;
        rows(r).Lost_Int = lostI;
        rows(r).Gained_Int = gainedI;

        rows(r).Base_Grp = bG;
        rows(r).Sensitivity_Grp = sG;
        rows(r).Lost_Grp = lostG;
        rows(r).Gained_Grp = gainedG;

        rows(r).Base_Cond = bC;
        rows(r).Sensitivity_Cond = sC;
        rows(r).Lost_Cond = lostC;
        rows(r).Gained_Cond = gainedC;

        rows(r).Any_FDR_Change = anyFDRchange;
        rows(r).Direction_Changed = dirChanged;

        rows(r).Pairwise_Family = string(famPW);
        rows(r).Pairwise_N_Matched = nMatchPW;
        rows(r).Pairwise_Jaccard = jacc;
        rows(r).Effect_Spearman = rho;
        rows(r).Effect_MedianAbsDiff = mad;

        rows(r).EMM_Family = string(famEMM);
        rows(r).EMM_N_Matched = nMatchEMM;
        rows(r).EMM_RMSE_Transformed = rmseT;
        rows(r).EMM_RMSE_BackTransformed = rmseBT;

    end
end

detail = struct2table(rows);

summary = summarize_extended(detail, ERPs, alpha);

OUT = struct();
OUT.detail = detail;
OUT.summary = summary;

end

% =============================== SUMMARY ===================================

function summary = summarize_extended(detail, ERPs, alpha)

% Overall + per-ERP rows with:
%   - counts (lost/gained/direction change)
%   - agreement + kappa for omnibus tests
%   - aggregated pairwise & EMM metrics

rows = table();

rows = [rows; build_summary_row(detail, "Overall", alpha)];

for e = 1: numel(ERPs)
    erpName = string(ERPs{e});
    sub = detail(detail.ERP == erpName, :);
    rows = [rows; build_summary_row(sub, erpName, alpha)];
end

summary = rows;

end

function t = build_summary_row(D, scopeName, alpha)

t = table();

t.Scope = string(scopeName);
t.Total_Models = height(D);

t.Any_FDR_Change = sum(D.Any_FDR_Change);
t.Direction_Changed = sum(D.Direction_Changed);

t.Lost_Any = sum(D.Lost_Int | D.Lost_Grp | D.Lost_Cond);
t.Gained_Any = sum(D.Gained_Int | D.Gained_Grp | D.Gained_Cond);

t.Lost_Int = sum(D.Lost_Int);
t.Gained_Int = sum(D.Gained_Int);

t.Lost_Grp = sum(D.Lost_Grp);
t.Gained_Grp = sum(D.Gained_Grp);

t.Lost_Cond = sum(D.Lost_Cond);
t.Gained_Cond = sum(D.Gained_Cond);

% ---------------- omnibus agreement + kappa ----------------

[bAgreeI, kI, nI, bPosI, sPosI] = agreement_kappa(D.Base_Int, D.Sensitivity_Int);
[bAgreeG, kG, nG, bPosG, sPosG] = agreement_kappa(D.Base_Grp, D.Sensitivity_Grp);
[bAgreeC, kC, nC, bPosC, sPosC] = agreement_kappa(D.Base_Cond, D.Sensitivity_Cond);

t.Omnibus_Int_N = nI;
t.Omnibus_Int_BaseSigN = bPosI;
t.Omnibus_Int_SensSigN = sPosI;
t.Omnibus_Int_AgreePct = bAgreeI;
t.Omnibus_Int_Kappa = kI;

t.Omnibus_Grp_N = nG;
t.Omnibus_Grp_BaseSigN = bPosG;
t.Omnibus_Grp_SensSigN = sPosG;
t.Omnibus_Grp_AgreePct = bAgreeG;
t.Omnibus_Grp_Kappa = kG;

t.Omnibus_Cond_N = nC;
t.Omnibus_Cond_BaseSigN = bPosC;
t.Omnibus_Cond_SensSigN = sPosC;
t.Omnibus_Cond_AgreePct = bAgreeC;
t.Omnibus_Cond_Kappa = kC;

% ---------------- pairwise metrics aggregates ----------------

t.Pairwise_CellsWithMetrics = sum(isfinite(D.Pairwise_Jaccard));

t.Pairwise_Jaccard_Median = median(D.Pairwise_Jaccard, 'omitnan');
t.Pairwise_Jaccard_Mean = mean(D.Pairwise_Jaccard, 'omitnan');

t.EffectRho_Median = median(D.Effect_Spearman, 'omitnan');
t.EffectRho_Mean = mean(D.Effect_Spearman, 'omitnan');

t.EffectMAD_Median = median(D.Effect_MedianAbsDiff, 'omitnan');
t.EffectMAD_Mean = mean(D.Effect_MedianAbsDiff, 'omitnan');

% ---------------- EMM RMSE aggregates ----------------

t.EMM_CellsWithMetrics = sum(isfinite(D.EMM_RMSE_Transformed));

t.EMM_RMSE_Trans_Median = median(D.EMM_RMSE_Transformed, 'omitnan');
t.EMM_RMSE_Trans_Mean = mean(D.EMM_RMSE_Transformed, 'omitnan');

t.EMM_RMSE_BT_Median = median(D.EMM_RMSE_BackTransformed, 'omitnan');
t.EMM_RMSE_BT_Mean = mean(D.EMM_RMSE_BackTransformed, 'omitnan');

t.AlphaFDR = alpha;

end

function [agreePct, kappa, n, nBasePos, nSensPos] = agreement_kappa(baseFlag, sensFlag)

baseFlag = logical(baseFlag(:));
sensFlag = logical(sensFlag(:));

n = numel(baseFlag);

nBasePos = sum(baseFlag);
nSensPos = sum(sensFlag);

agree = mean(baseFlag == sensFlag);
agreePct = 100 * agree;

pA1 = mean(baseFlag);
pB1 = mean(sensFlag);
pA0 = 1 - pA1;
pB0 = 1 - pB1;

pe = pA1 * pB1 + pA0 * pB0;

den = (1 - pe);

if den <= 0
    kappa = NaN;
else
    kappa = (agree - pe) / den;
end

end

% =============================== DETAIL METRICS ============================

function [jacc, rho, mad, fam, nMatch] = pairwise_concordance(Rb, Rs, alpha)

% Returns NaN when insufficient data.
% Correlation + MAD computed across ALL matched contrasts (not only significant).
% Jaccard computed on significant sets only (pFDR < alpha).

jacc = NaN;
rho = NaN;
mad = NaN;
fam = "";
nMatch = 0;

if ~isstruct(Rb) || ~isstruct(Rs)
    return;
end

bI = lt_alpha(safeget(Rb, NaN, 'tests', 'interaction', 'pFDR'), alpha);
sI = lt_alpha(safeget(Rs, NaN, 'tests', 'interaction', 'pFDR'), alpha);

if bI && sI
    Tb = collect_pairwise_interaction(Rb);
    Ts = collect_pairwise_interaction(Rs);
    fam = "interaction";
    [jacc, rho, mad, nMatch] = compare_pairwise_tables(Tb, Ts, alpha);
    return;
end

if ~bI && ~sI

    bG = lt_alpha(safeget(Rb, NaN, 'tests', 'group', 'pFDR'), alpha);
    sG = lt_alpha(safeget(Rs, NaN, 'tests', 'group', 'pFDR'), alpha);

    if bG && sG
        Tb = collect_pairwise_main(Rb, 'Group');
        Ts = collect_pairwise_main(Rs, 'Group');
        fam = "group";
        [jacc, rho, mad, nMatch] = compare_pairwise_tables(Tb, Ts, alpha);
        return;
    end

    bC = lt_alpha(safeget(Rb, NaN, 'tests', 'condition', 'pFDR'), alpha);
    sC = lt_alpha(safeget(Rs, NaN, 'tests', 'condition', 'pFDR'), alpha);

    if bC && sC
        Tb = collect_pairwise_main(Rb, 'Condition');
        Ts = collect_pairwise_main(Rs, 'Condition');
        fam = "condition";
        [jacc, rho, mad, nMatch] = compare_pairwise_tables(Tb, Ts, alpha);
        return;
    end

end

end

function [rmseT, rmseBT, fam, nMatch] = emm_rmse_concordance(Rb, Rs)

% RMSE between EMM vectors across pipelines, within the family actually reported.
% Computed only when BOTH pipelines have relevant EMM tables.

rmseT = NaN;
rmseBT = NaN;
fam = "";
nMatch = 0;

if ~isstruct(Rb) || ~isstruct(Rs)
    return;
end

pB_I = safeget(Rb, NaN, 'tests', 'interaction', 'pFDR');
pS_I = safeget(Rs, NaN, 'tests', 'interaction', 'pFDR');

bI = isfinite(pB_I) && pB_I < 0.05;
sI = isfinite(pS_I) && pS_I < 0.05;

if bI && sI
    EB = get_first_table(safeget(Rb, struct(), 'emm'), {'Group_by_Condition', 'Condition_by_Group'});
    ES = get_first_table(safeget(Rs, struct(), 'emm'), {'Group_by_Condition', 'Condition_by_Group'});

    if isempty(EB) || isempty(ES)
        return;
    end

    fam = "interaction";

    [rmseT, rmseBT, nMatch] = rmse_emm_tables(EB, ES, {'Group', 'Condition'});
    return;
end

if ~bI && ~sI

    EB = [];
    ES = [];

    pB_G = safeget(Rb, NaN, 'tests', 'group', 'pFDR');
    pS_G = safeget(Rs, NaN, 'tests', 'group', 'pFDR');

    if isfinite(pB_G) && isfinite(pS_G) && pB_G < 0.05 && pS_G < 0.05
        EB = safeget(Rb, [], 'emm', 'Group');
        ES = safeget(Rs, [], 'emm', 'Group');
        fam = "group";
        if ~isempty(EB) && ~isempty(ES)
            [rmseT, rmseBT, nMatch] = rmse_emm_tables(EB, ES, {'Group'});
            return;
        end
    end

    pB_C = safeget(Rb, NaN, 'tests', 'condition', 'pFDR');
    pS_C = safeget(Rs, NaN, 'tests', 'condition', 'pFDR');

    if isfinite(pB_C) && isfinite(pS_C) && pB_C < 0.05 && pS_C < 0.05
        EB = safeget(Rb, [], 'emm', 'Condition');
        ES = safeget(Rs, [], 'emm', 'Condition');
        fam = "condition";
        if ~isempty(EB) && ~isempty(ES)
            [rmseT, rmseBT, nMatch] = rmse_emm_tables(EB, ES, {'Condition'});
            return;
        end
    end

end

end

function [rmseT, rmseBT, nMatch] = rmse_emm_tables(Tb, Ts, keyVars)

rmseT = NaN;
rmseBT = NaN;
nMatch = 0;

if ~istable(Tb) || ~istable(Ts)
    return;
end

Kb = emm_keys(Tb, keyVars);
Ks = emm_keys(Ts, keyVars);

if isempty(Kb) || isempty(Ks)
    return;
end

[common, ib, is] = intersect(Kb, Ks, 'stable');

if numel(common) < 2
    return;
end

xb = pick_col(Tb, {'Estimate'});
xs = pick_col(Ts, {'Estimate'});

db = xb(ib);
ds = xs(is);

good = isfinite(db) & isfinite(ds);

if sum(good) >= 2
    nMatch = sum(good);
    rmseT = sqrt(mean((db(good) - ds(good)).^2, 'omitnan'));
end

% Back-transformed (optional)
if ismember('Estimate_BT', Tb.Properties.VariableNames) && ismember('Estimate_BT', Ts.Properties.VariableNames)

    bb = Tb.Estimate_BT(ib);
    ss = Ts.Estimate_BT(is);

    good2 = isfinite(bb) & isfinite(ss);

    if sum(good2) >= 2
        rmseBT = sqrt(mean((bb(good2) - ss(good2)).^2, 'omitnan'));
    end

end

end

function K = emm_keys(T, keyVars)

K = strings(height(T), 1);

for i = 1: height(T)

    parts = strings(1, numel(keyVars));

    for k = 1: numel(keyVars)

        vn = keyVars{k};

        if ~ismember(vn, T.Properties.VariableNames)
            parts(k) = "";
        else
            parts(k) = string(T.(vn)(i));
        end

    end

    K(i) = strjoin(parts, "|");

end

end

function [jacc, rho, mad, nMatch] = compare_pairwise_tables(Tb, Ts, alpha)

jacc = NaN;
rho = NaN;
mad = NaN;
nMatch = 0;

if isempty(Tb) || isempty(Ts) || ~istable(Tb) || ~istable(Ts)
    return;
end

if ~ismember('pFDR', Tb.Properties.VariableNames) || ~ismember('pFDR', Ts.Properties.VariableNames)
    return;
end

Kb = pairwise_keys(Tb);
Ks = pairwise_keys(Ts);

if isempty(Kb) || isempty(Ks)
    return;
end

sigB = isfinite(Tb.pFDR) & (Tb.pFDR < alpha);
sigS = isfinite(Ts.pFDR) & (Ts.pFDR < alpha);

setB = unique(Kb(sigB));
setS = unique(Ks(sigS));

if isempty(setB) && isempty(setS)
    jacc = 1;
else
    inter = intersect(setB, setS);
    uni = union(setB, setS);
    jacc = numel(inter) / max(1, numel(uni));
end

if ~ismember('Estimate', Tb.Properties.VariableNames) || ~ismember('Estimate', Ts.Properties.VariableNames)
    return;
end

[common, ib, is] = intersect(Kb, Ks, 'stable');

if numel(common) < 3
    return;
end

eb = Tb.Estimate(ib);
es = Ts.Estimate(is);

good = isfinite(eb) & isfinite(es);

if sum(good) < 3
    return;
end

nMatch = sum(good);

rho = spearman_robust(eb(good), es(good));
mad = median(abs(eb(good) - es(good)), 'omitnan');

end

function rho = spearman_robust(x, y)

rho = NaN;

try
    rho = corr(x(:), y(:), 'Type', 'Spearman');
catch
    try
        rx = tiedrank(x(:));
        ry = tiedrank(y(:));
        C = corrcoef(rx, ry);
        rho = C(1, 2);
    catch
        rho = NaN;
    end
end

end

function T = collect_pairwise_interaction(R)

T = table();

if ~isfield(R, 'pairwise') || isempty(R.pairwise)
    return;
end

parts = {};

P1 = get_first_table(R.pairwise, {'Group_within_Condition', 'Group_by_Condition'});
P2 = get_first_table(R.pairwise, {'Condition_within_Group', 'Condition_by_Group'});

if ~isempty(P1)
    P1.Family = repmat("Group_within_Condition", height(P1), 1);
    parts{end + 1} = P1;
end

if ~isempty(P2)
    P2.Family = repmat("Condition_within_Group", height(P2), 1);
    parts{end + 1} = P2;
end

if ~isempty(parts)
    T = vertcat(parts{:});
end

end

function T = collect_pairwise_main(R, which)

T = table();

if ~isfield(R, 'pairwise') || isempty(R.pairwise)
    return;
end

if strcmpi(which, 'Group')

    if isfield(R.pairwise, 'Group') && istable(R.pairwise.Group)
        T = R.pairwise.Group;
        T.Family = repmat("Group_main", height(T), 1);
    end

elseif strcmpi(which, 'Condition')

    if isfield(R.pairwise, 'Condition') && istable(R.pairwise.Condition)
        T = R.pairwise.Condition;
        T.Family = repmat("Condition_main", height(T), 1);
    end

end

end

function K = pairwise_keys(P)

% Key = Family | context | contrast
%   - context: Condition=... (if present) else Group=... (if present) else ""
%   - contrast: prefer the contrast column (Group or Condition) if present

K = strings(height(P), 1);

fam = repmat("", height(P), 1);

if ismember('Family', P.Properties.VariableNames)
    fam = string(P.Family);
end

ctx = repmat("", height(P), 1);

if ismember('Condition', P.Properties.VariableNames)
    ctx = "Condition=" + string(P.Condition);
elseif ismember('Group', P.Properties.VariableNames)
    ctx = "Group=" + string(P.Group);
end

contrast = repmat("", height(P), 1);

% Prefer explicit contrast column:
%   - In your lme_posthoc tables, the contrast label is usually stored in a column
%     named exactly as the factor being contrasted ('Group' or 'Condition').
if ismember('Group', P.Properties.VariableNames) && ismember('Condition', P.Properties.VariableNames)
    % For interaction tables, the contrasted factor is whichever column contains "A - B".
    % Heuristic: pick the column that contains " - " in at least one row.
    gStr = string(P.Group);
    cStr = string(P.Condition);
    if any(contains(gStr, " - "))
        contrast = gStr;
    elseif any(contains(cStr, " - "))
        contrast = cStr;
    else
        contrast = gStr;
    end
elseif ismember('Group', P.Properties.VariableNames)
    contrast = string(P.Group);
elseif ismember('Condition', P.Properties.VariableNames)
    contrast = string(P.Condition);
elseif ismember('Contrast', P.Properties.VariableNames)
    contrast = string(P.Contrast);
else
    % fallback: first string/categorical column
    for j = 1: width(P)
        if isstring(P{:, j}) || iscategorical(P{:, j})
            contrast = string(P{:, j});
            break;
        end
    end
end

for i = 1: height(P)
    K(i) = fam(i) + "|" + ctx(i) + "|" + contrast(i);
end

end

% =============================== UTILITIES =================================

function alpha = pick_alpha(Out_base, Out_sens)

alpha = 0.05;

a1 = safeget(Out_base, NaN, 'info', 'alphaFDR');
a2 = safeget(Out_sens, NaN, 'info', 'alphaFDR');

if isfinite(a1)
    alpha = a1;
elseif isfinite(a2)
    alpha = a2;
end

end

function measureField = resolve_measure_field(Out_base, Out_sens, measure)

if strcmpi(measure, 'Occurrence')

    if isfield(Out_base, 'Num_occurrence') && isfield(Out_sens, 'Num_occurrence')
        measureField = 'Num_occurrence';
        return;
    end

    if isfield(Out_base, 'Occurrence') && isfield(Out_sens, 'Occurrence')
        measureField = 'Occurrence';
        return;
    end

    if isfield(Out_base, 'Num_occurrence')
        measureField = 'Num_occurrence';
        return;
    end

    measureField = 'Occurrence';
    return;
end

measureField = char(measure);

end

function tf = lt_alpha(p, alpha)

if ~isfinite(p)
    tf = false;
else
    tf = (p < alpha);
end

end

function val = safeget(S, defaultVal, varargin)

val = defaultVal;

try
    for k = 1: numel(varargin)
        S = S.(varargin{k});
    end
    val = S;
catch
end

end

function T = get_first_table(S, candidates)

T = [];

if isempty(S) || ~isstruct(S)
    return;
end

for k = 1: numel(candidates)
    fn = candidates{k};
    if isfield(S, fn) && istable(S.(fn)) && ~isempty(S.(fn))
        T = S.(fn);
        return;
    end
end

end

function x = pick_col(T, candidates)

x = NaN(height(T), 1);

for k = 1: numel(candidates)
    cn = candidates{k};
    if ismember(cn, T.Properties.VariableNames)
        x = T.(cn);
        return;
    end
end

end

function ok = same_ordering_1F(T1, T2, factorName)

ok = true;

if isempty(T1) || isempty(T2) || ~istable(T1) || ~istable(T2)
    return;
end

if ~ismember(factorName, string(T1.Properties.VariableNames)) || ~ismember(factorName, string(T2.Properties.VariableNames))
    return;
end

y1 = pick_emm_values(T1);
y2 = pick_emm_values(T2);

lv1 = string(T1.(factorName));
lv2 = string(T2.(factorName));

common = intersect(lv1, lv2, 'stable');

if numel(common) < 2
    return;
end

a = nan(numel(common), 1);
b = nan(numel(common), 1);

for i = 1: numel(common)
    a(i) = y1(find(lv1 == common(i), 1, 'first'));
    b(i) = y2(find(lv2 == common(i), 1, 'first'));
end

[~, ordA] = sort(a, 'ascend');
[~, ordB] = sort(b, 'ascend');

ok = isequal(ordA, ordB);

end

function ok = same_ordering_2F(T1, T2, factorA, factorB)

ok = true;

if isempty(T1) || isempty(T2) || ~istable(T1) || ~istable(T2)
    return;
end

need1 = ismember(string({factorA, factorB}), string(T1.Properties.VariableNames));
need2 = ismember(string({factorA, factorB}), string(T2.Properties.VariableNames));

if ~all(need1) || ~all(need2)
    return;
end

bLv = intersect(string(unique(T1.(factorB))), string(unique(T2.(factorB))), 'stable');

if isempty(bLv)
    return;
end

for k = 1: numel(bLv)

    bk = bLv(k);

    sub1 = T1(string(T1.(factorB)) == bk, :);
    sub2 = T2(string(T2.(factorB)) == bk, :);

    if height(sub1) < 2 || height(sub2) < 2
        continue;
    end

    ysub1 = pick_emm_values(sub1);
    ysub2 = pick_emm_values(sub2);

    lv1 = string(sub1.(factorA));
    lv2 = string(sub2.(factorA));

    common = intersect(lv1, lv2, 'stable');

    if numel(common) < 2
        continue;
    end

    a = nan(numel(common), 1);
    b = nan(numel(common), 1);

    for i = 1: numel(common)
        a(i) = ysub1(find(lv1 == common(i), 1, 'first'));
        b(i) = ysub2(find(lv2 == common(i), 1, 'first'));
    end

    [~, ordA] = sort(a, 'ascend');
    [~, ordB] = sort(b, 'ascend');

    if ~isequal(ordA, ordB)
        ok = false;
        return;
    end

end

end

function y = pick_emm_values(T)

if ismember('Estimate', T.Properties.VariableNames)
    y = T.Estimate;
elseif ismember('Estimate_BT', T.Properties.VariableNames)
    y = T.Estimate_BT;
elseif ismember('EMM_BT', T.Properties.VariableNames)
    y = T.EMM_BT;
else
    y = nan(height(T), 1);
end

end
