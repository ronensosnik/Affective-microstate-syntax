function OUT = lme_posthoc(ALL_Mean, measure)

% Fits LME models, computes EMMs, and pairwise contrasts (direction & CI)
% without using margins(). EMMs are back-transformed to ms / % / Hz.

measure = validatestring(measure, {'Duration', 'Coverage', 'Num_occurrence'});

ERPs = {'N200', 'P300', 'LPP'};
Micros = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};

winSec = struct('N200', 0.12, 'P300', 0.2, 'LPP', 0.5);
winSamp = struct('N200', 120, 'P300', 200, 'LPP', 500);

OUT = struct();
OUT.(measure) = cell(3, 7);

p_int = nan(21, 1);
p_grp = nan(21, 1);
p_cond = nan(21, 1);
map = nan(21, 2);
row = 0;

for e = 1: 3

    ERPname = ERPs{e};
    wS = winSec.(ERPname);
    wN = winSamp.(ERPname);

    epsC = 0.5 / wN; % clamp for coverage fraction (avoid 0/1)
    cOcc = 0.5 / wS; % continuity correction in Hz (0.5 / window-seconds)

    for m = 1: 7

        row = row + 1;
        map(row, :) = [e m];
        S = ALL_Mean.(ERPname).(Micros{m});

        % -------- Build analysis table --------

        subj = categorical(S.Number(:));

        gStr = string(S.Group(:));
        gU = unique(gStr, 'stable');

        if any(gU == "HC")
            gOrder = ["HC"; setdiff(gU, "HC", 'stable')];
        else
            gOrder = gU;
        end

        T = table();
        T.SubjectID = subj;
        T.Group = categorical(gStr, gOrder);
        T.Condition = categorical(string(S.Condition(:)), ["Neutral" "Positive" "Negative"]);
        T.Age = double(S.Age(:));

        % -------- Transform response (static metrics) --------

        switch measure

            case 'Duration'
                resp = double(S.Duration(:)); % ms
                trans = struct('type', 'none', 'respName', 'Duration', 'ERP', ERPname, 'Microstate', Micros{m});

            case 'Coverage'
                % Dataset Coverage is in percent: convert to fraction for logit.
                p = double(S.Coverage(:)) / 100;

                % Window-specific clamp at 0.5 / number-of-samples (and symmetric upper clamp).
                p = min(max(p, epsC), 1 - epsC);

                resp = log(p ./ (1 - p)); % logit(fraction)
                trans = struct('type', 'logit', 'respName', 'Coverage_logit', 'ERP', ERPname, 'Microstate', Micros{m}, 'epsC', epsC, 'winSamp', wN);

            case 'Num_occurrence'
                % Occurrence is a non-negative rate in Hz: apply additive continuity correction then log.
                x = double(S.Num_occurrence(:));
                xCorr = x + cOcc;

                resp = log(xCorr);
                trans = struct('type', 'log', 'respName', 'Occurrence_log', 'ERP', ERPname, 'Microstate', Micros{m}, 'cOcc', cOcc, 'winSec', wS);

        end

        T.Response = resp;

        % Drop missing / undefined rows explicitly for stability
        keep = isfinite(T.Response) & isfinite(T.Age) & ~isundefined(T.Group) & ~isundefined(T.Condition) & ~isundefined(T.SubjectID);
        T = T(keep, :);

        % Center age within this ERP × microstate analysis
        T.Age_c = T.Age - mean(T.Age, 'omitnan');

        % -------- Random-effects structure: try random slope for Condition --------

        reSlope = '(1 + Condition | SubjectID)';
        reInt = '(1 | SubjectID)';

        useSlope = true;
        fTest = sprintf('Response ~ Group * Condition + Age_c + %s', reSlope);

        try
            m_tmp = fitlme(T, fTest, 'FitMethod', 'ML');
            if is_singular_cov(m_tmp)
                useSlope = false;
            end
        catch
            useSlope = false;
        end

        if useSlope
            reStr = reSlope;
        else
            reStr = reInt;
        end

        % -------- Omnibus models (ML for LRT) --------

        f_full = sprintf('Response ~ Group * Condition + Age_c + %s', reStr);
        f_no_int = sprintf('Response ~ Group + Condition + Age_c + %s', reStr);
        f_no_grp = sprintf('Response ~ Condition + Age_c + %s', reStr);
        f_no_cond = sprintf('Response ~ Group + Age_c + %s', reStr);

        m_full_int = fitlme(T, f_full, 'FitMethod', 'ML');
        m_no_int = fitlme(T, f_no_int, 'FitMethod', 'ML');
        m_no_grp = fitlme(T, f_no_grp, 'FitMethod', 'ML');
        m_no_cond = fitlme(T, f_no_cond, 'FitMethod', 'ML');

        [p_int(row), ~, ~] = safe_lrt_p(m_no_int, m_full_int);
        [p_grp(row), ~, ~] = safe_lrt_p(m_no_grp, m_no_int);
        [p_cond(row), ~, ~] = safe_lrt_p(m_no_cond, m_no_int);

        % Store raw objects

        R = struct();
        R.dataTable = T;
        R.transform = trans;
        R.randomEffects = struct('usedRandomSlopeForCondition', useSlope, 'reStr', reStr);

        R.models.full_with_interaction = m_full_int;
        R.models.no_interaction_ML = m_no_int;
        R.models.no_group = m_no_grp;
        R.models.no_condition = m_no_cond;

        R.tests = struct();
        R.tests.interaction = struct('pRaw', p_int(row));
        R.tests.group = struct('pRaw', p_grp(row));
        R.tests.condition = struct('pRaw', p_cond(row));

        R.emm = struct();
        R.pairwise = struct();
        R.model_used = '';

        OUT.(measure){e, m} = R;

    end
end

% -------- Per-ERP BH-FDR across omnibus tests (7 microstates × 3 tests per ERP) --------

q_int = nan(21, 1);
q_grp = nan(21, 1);
q_cond = nan(21, 1);

for ee = 1: 3
    idx = find(map(:, 1) == ee);
    n = numel(idx);

    pStack = [p_int(idx); p_grp(idx); p_cond(idx)];
    keep = isfinite(pStack);

    qStack = nan(size(pStack));

    if any(keep)
        qStack(keep) = bh_adjust(pStack(keep));
    end

    q_int(idx) = qStack(1: n);
    q_grp(idx) = qStack(n + 1: 2 * n);
    q_cond(idx) = qStack(2 * n + 1: 3 * n);
end

% -------- EMMs and Pairwise contrasts (no margins dependency) --------

for r = 1: 21
    e = map(r, 1);
    m = map(r, 2);

    R = OUT.(measure){e, m};
    T = R.dataTable;
    reStr = R.randomEffects.reStr;

    R.tests.interaction.pRaw = p_int(r);
    R.tests.group.pRaw = p_grp(r);
    R.tests.condition.pRaw = p_cond(r);

    R.tests.interaction.pFDR = q_int(r);
    R.tests.group.pFDR = q_grp(r);
    R.tests.condition.pFDR = q_cond(r);

    atSpec = struct('Age_c', 0);

    if isfinite(q_int(r)) && q_int(r) < 0.05

        % ----- Interaction present: simple-effects EMMs + pairwise -----

        R.model_used = 'full_with_interaction';

        R.emm.Group_by_Condition = emm_lme(R.models.full_with_interaction, {'Group', 'Condition'}, T, atSpec);
        R.emm.Condition_by_Group = emm_lme(R.models.full_with_interaction, {'Condition', 'Group'}, T, atSpec);

        R.emm.Group_by_Condition = addBT_emm(R.emm.Group_by_Condition, measure, R.transform);
        R.emm.Condition_by_Group = addBT_emm(R.emm.Condition_by_Group, measure, R.transform);

        f_full_ref = sprintf('Response ~ Group * Condition + Age_c + %s', reStr);
        m_full_int_ref = fitlme(T, f_full_ref, 'DummyVarCoding', 'reference', 'FitMethod', 'ML');

        [PG, famG] = pairwise_simple(m_full_int_ref, T, measure, 'Group', 'Condition');
        PG.pFDR = within_families_fdr(PG.pValue, famG);
        R.pairwise.Group_within_Condition = PG;

        [PC, famC] = pairwise_simple(m_full_int_ref, T, measure, 'Condition', 'Group');
        PC.pFDR = within_families_fdr(PC.pValue, famC);
        R.pairwise.Condition_within_Group = PC;

    else

        % ----- No interaction: main-effects EMMs + pairwise -----

        R.model_used = 'no_interaction';

        f_main_ref = sprintf('Response ~ Group + Condition + Age_c + %s', reStr);
        m_main_ref = fitlme(T, f_main_ref, 'DummyVarCoding', 'reference', 'FitMethod', 'REML');

        if isfinite(q_grp(r)) && q_grp(r) < 0.05

            R.emm.Group = emm_lme(m_main_ref, {'Group'}, T, atSpec);
            R.emm.Group = addBT_emm(R.emm.Group, measure, R.transform);

            PG = pairwise_main(m_main_ref, T, measure, 'Group');
            PG.pFDR = bh_adjust(PG.pValue);

            R.pairwise.Group = PG;

        end

        if isfinite(q_cond(r)) && q_cond(r) < 0.05

            R.emm.Condition = emm_lme(m_main_ref, {'Condition'}, T, atSpec);
            R.emm.Condition = addBT_emm(R.emm.Condition, measure, R.transform);

            PC = pairwise_main(m_main_ref, T, measure, 'Condition');
            PC.pFDR = bh_adjust(PC.pValue);

            R.pairwise.Condition = PC;

        end

        R.models.no_interaction_REML = m_main_ref;

    end

    OUT.(measure){e, m} = R;
end

% -------- Global info --------

OUT.info.ERPs = ERPs;
OUT.info.Microstates = Micros;
OUT.info.measure = measure;
OUT.info.alphaFDR = 0.05;
OUT.info.correction = 'BH-FDR per ERP across 7 microstates × 3 omnibus tests (interaction, group, condition) for this measure';
OUT.info.notes = 'Coverage modeled as logit(fraction) with clamp 0.5 / window-samples; occurrence modeled as log(rate + 0.5 / window-seconds); duration modeled on identity. EMMs predicted with Conditional = false and back-transformed. Random slope for Condition attempted and dropped if failing/singular.';

end % ================= END MAIN =================


function tf = is_singular_cov(M)

% Heuristic singularity check for random-effects covariance parameters.
% Flags near-zero variance components, which often indicate singular fits.

tf = false;

try
    cp = M.CovarianceParameters;
    if istable(cp) && ismember('Estimate', cp.Properties.VariableNames)
        est = cp.Estimate;
        if any(~isfinite(est))
            tf = true;
            return;
        end
        if any(est < 1e-6)
            tf = true;
            return;
        end
    end
catch
    tf = false;
end

end


function tbl = emm_lme(M, factors, Tfit, at)

% Build a grid over "factors" and predict fixed-effects means (CIs).

if nargin < 4 || isempty(at)
    at = struct;
end

nF = numel(factors);
levels = cell(nF, 1);
sz = zeros(1, nF);

for i = 1: nF
    levels{i} = categories(Tfit.(factors{i}));
    sz(i) = numel(levels{i});
end

subs = cell(1, nF);
[subs{:}] = ndgrid(levels{:});
N = numel(subs{1});

newT = table;

for i = 1: nF
    lv = string(subs{i}(:));
    newT.(factors{i}) = categorical(lv, categories(Tfit.(factors{i})));
end

newT.SubjectID = repmat(Tfit.SubjectID(1), N, 1);

vars = setdiff(M.PredictorNames, [factors, {'SubjectID'}]);

for i = 1: numel(vars)
    vn = vars{i};

    if ismember(vn, newT.Properties.VariableNames)
        continue;
    end

    if isfield(at, vn)
        newT.(vn) = repmat(at.(vn), N, 1);
    else
        if isnumeric(Tfit.(vn))
            newT.(vn) = repmat(mean(Tfit.(vn), 'omitnan'), N, 1);
        else
            newT.(vn) = repmat(Tfit.(vn)(1), N, 1);
        end
    end
end

[yhat, yCI] = predict(M, newT, 'Conditional', false);

tbl = newT(:, factors);
tbl.Estimate = yhat;
tbl.CI_lo = yCI(:, 1);
tbl.CI_hi = yCI(:, 2);

end


function tbl = addBT_emm(tbl, measure, trans)

switch trans.type

    case 'none'
        tbl.Estimate_BT = tbl.Estimate;
        tbl.CI_lo_BT = tbl.CI_lo;
        tbl.CI_hi_BT = tbl.CI_hi;
        tbl.Units = repmat("ms", height(tbl), 1);

    case 'logit'
        % Back-transform to percent scale.
        tbl.Estimate_BT = 100 * (1 ./ (1 + exp(-tbl.Estimate)));
        tbl.CI_lo_BT = 100 * (1 ./ (1 + exp(-tbl.CI_lo)));
        tbl.CI_hi_BT = 100 * (1 ./ (1 + exp(-tbl.CI_hi)));
        tbl.Units = repmat("%", height(tbl), 1);

    case 'log'
        % Back-transform to rate in Hz, removing the continuity correction.
        cOcc = trans.cOcc;

        tbl.Estimate_BT = exp(tbl.Estimate) - cOcc;
        tbl.CI_lo_BT = exp(tbl.CI_lo) - cOcc;
        tbl.CI_hi_BT = exp(tbl.CI_hi) - cOcc;

        % Numerical safety: rates cannot be negative.
        tbl.Estimate_BT = max(tbl.Estimate_BT, 0);
        tbl.CI_lo_BT = max(tbl.CI_lo_BT, 0);
        tbl.CI_hi_BT = max(tbl.CI_hi_BT, 0);

        tbl.Units = repmat("Hz", height(tbl), 1);

    otherwise
        tbl.Estimate_BT = tbl.Estimate;
        tbl.CI_lo_BT = tbl.CI_lo;
        tbl.CI_hi_BT = tbl.CI_hi;
        tbl.Units = repmat("", height(tbl), 1);

end

tbl.Measure = repmat(string(measure), height(tbl), 1);

end


function [tbl, fam] = pairwise_simple(Mref, T, measure, A, B)

lvA = categories(T.(A));
lvB = categories(T.(B));

pairs = nchoosek(1: numel(lvA), 2);

rows = [];

for jb = 1: numel(lvB)

    bLevel = lvB{jb};

    for i = 1: size(pairs, 1)

        a1 = lvA{pairs(i, 1)};
        a2 = lvA{pairs(i, 2)};

        c = fixedContrast_simple(Mref, T, A, B, a1, a2, bLevel, measure);
        rows = [rows; c];

    end
end

tbl = struct2table(rows);
fam = tbl.(B);

end


function c = fixedContrast_simple(Mref, T, A, B, a1, a2, bLevel, measure)

at = struct('Age_c', 0);
tbl1 = emm_lme(Mref, {A, B}, T, at);

idx1 = (string(tbl1.(A)) == string(a1)) & (string(tbl1.(B)) == string(bLevel));
idx2 = (string(tbl1.(A)) == string(a2)) & (string(tbl1.(B)) == string(bLevel));

b1 = tbl1.Estimate(idx1);
b2 = tbl1.Estimate(idx2);

se1 = (tbl1.CI_hi(idx1) - tbl1.Estimate(idx1)) / 1.96;
se2 = (tbl1.CI_hi(idx2) - tbl1.Estimate(idx2)) / 1.96;

diff = b1 - b2;
se = sqrt(se1.^2 + se2.^2);
z = diff ./ se;
p = 2 * (1 - normcdf(abs(z)));

c = struct();
c.(A) = string(a1) + " - " + string(a2);
c.(B) = string(bLevel);
c.Estimate = diff;
c.SE = se;
c.zValue = z;
c.pValue = p;
c.Measure = string(measure);

end


function tbl = pairwise_main(Mref, T, measure, A)

lv = categories(T.(A));
pairs = nchoosek(1: numel(lv), 2);

rows = [];

for i = 1: size(pairs, 1)

    a1 = lv{pairs(i, 1)};
    a2 = lv{pairs(i, 2)};

    c = fixedContrast_main(Mref, T, A, a1, a2);
    rows = [rows; c];

end

tbl = struct2table(rows);
tbl.Measure = repmat(string(measure), height(tbl), 1);

end


function c = fixedContrast_main(Mref, T, A, a1, a2)

at = struct('Age_c', 0);
tbl1 = emm_lme(Mref, {A}, T, at);

idx1 = (string(tbl1.(A)) == string(a1));
idx2 = (string(tbl1.(A)) == string(a2));

b1 = tbl1.Estimate(idx1);
b2 = tbl1.Estimate(idx2);

se1 = (tbl1.CI_hi(idx1) - tbl1.Estimate(idx1)) / 1.96;
se2 = (tbl1.CI_hi(idx2) - tbl1.Estimate(idx2)) / 1.96;

diff = b1 - b2;
se = sqrt(se1.^2 + se2.^2);
z = diff ./ se;
p = 2 * (1 - normcdf(abs(z)));

c = struct();
c.(A) = string(a1) + " - " + string(a2);
c.Estimate = diff;
c.SE = se;
c.zValue = z;
c.pValue = p;

end


function p_adj = within_families_fdr(p, fam)

p_adj = nan(size(p));
u = unique(fam);

for k = 1: numel(u)
    idx = (fam == u(k));
    p_adj(idx) = bh_adjust(p(idx));
end

end


function q = bh_adjust(p)

p = p(:);
keep = isfinite(p);

q = nan(size(p));

if ~any(keep)
    return;
end

pp = p(keep);
[ps, ord] = sort(pp, 'ascend');
m = numel(ps);

qtmp = ps .* m ./ (1: m)';
qtmp = flipud(cummin(flipud(qtmp)));

qtmp(qtmp > 1) = 1;

qkeep = nan(size(pp));
qkeep(ord) = qtmp;

q(keep) = qkeep;

end


function [p, LR, df, dLL, method] = safe_lrt_p(m_reduced, m_full)

% Robust likelihood-ratio test p-value for nested LinearMixedModel objects.

p = NaN;
LR = NaN;
df = NaN;
dLL = NaN;
method = "";

try
    cmp = compare(m_reduced, m_full, 'CheckNesting', true);
    p = cmp.pValue(2);
    LR = cmp.LRStat(2);
    df = cmp.DF(2);
    dLL = 0.5 * LR;
    method = "compare";
    return;
catch
end

try
    ll_full = m_full.LogLikelihood;
    ll_red = m_reduced.LogLikelihood;
    dLL = ll_full - ll_red;

    if ~isfinite(dLL)
        p = NaN;
        method = "manual_nan";
        return;
    end

    LR = 2 * max(0, dLL);

    df = size(m_full.Coefficients, 1) - size(m_reduced.Coefficients, 1);

    if ~isfinite(df) || df <= 0
        df = 1;
    end

    p = 1 - chi2cdf(LR, df);

    if ~isfinite(p)
        p = NaN;
        method = "manual_nan";
    else
        p = max(0, min(1, p));
        method = "manual";
    end
catch
    p = NaN;
    method = "failed";
end

end
