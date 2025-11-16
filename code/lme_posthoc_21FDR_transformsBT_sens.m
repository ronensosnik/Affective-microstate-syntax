function OUT = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, measure, alphaFDR, medSpec)
% LME_POSTHOC_21FDR_TRANSFORMSBT (revised: adds medication covariates)
% Fits LME models, computes EMMs, and pairwise contrasts (direction & CI)
% without using margins(). EMMs are back-transformed to ms / % / Hz.
%
% Usage (primary medication adjustment with five class dummies):
%   OUT = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05);
%
% Usage (sensitivity, low DF: poly_count only):
%   OUT = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05, 'poly');

if nargin < 3
    alphaFDR = 0.05;
end

if nargin < 4 || isempty(medSpec)
    medSpec = 'classes'; % 'classes' (AD,AP,MS,ANX,OTHER) or 'poly' (poly_count only)
end

medSpec = validatestring(medSpec, {'classes', 'poly'});
measure = validatestring(measure, {'Duration', 'Coverage', 'Occurrence'});

ERPs = {'N200', 'P300', 'LPP'};
Micros = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};

switch measure
    case 'Duration'
        respField = 'Duration';
        respName= 'Duration';

    case 'Coverage'
        respField = 'Coverage';
        respName = 'logitCov';

    case 'Occurrence'
        respField = 'Num_occurrence';
        respName = 'logOcc';
end

winSec = struct('N200', 0.120, 'P300', 0.200, 'LPP', 0.500);
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
    cOcc = 0.5 / wS; % continuity for occurrence (Hz)

    for m = 1: 7

        row = row + 1;
        map(row,:) = [e m];
        S = ALL_Mean.(ERPname).(Micros{m});

        % -------- Build analysis table & transform response --------

        T = table(S.Number(:), categorical(S.Group(:)), categorical(S.Condition(:)), S.Age(:), 'VariableNames', {'SubjectID', 'Group', 'Condition', 'Age'});

        switch measure
            case 'Duration'
                resp = S.Duration(:); % ms

            case 'Coverage'
                p = S.Coverage(:) / 100;
                p = min(max(p,epsC),1 - epsC);
                resp = log(p./(1 - p)); % logit

            case 'Occurrence'
                resp = log(S.Num_occurrence(:) + cOcc); % log(rate + c)
        end

        T.(respName) = resp;
        T.SubjectID = categorical(T.SubjectID);
        T.Age_c = T.Age - mean(T.Age, 'omitnan');

        % -------- Medications: add covariates from S.Meds (cell of structs) --------

        [MedTbl, medNames] = meds_to_table(S.Meds, medSpec);
        T = [T, MedTbl];

        T = rmmissing(T, 'DataVariables', {respName});

        % -------- Omnibus models (effects coding for valid LRTs) --------

        opts = {'FitMethod', 'ML', 'DummyVarCoding', 'effects'};
        medTerm = strjoin(medNames, ' + ');

        try
            fm_full_int = sprintf('%s ~ Group*Condition + Age_c + %s + (1 + Condition | SubjectID)', respName, medTerm);
            fm_full_noint = sprintf('%s ~ Group + Condition + Age_c + %s + (1 + Condition | SubjectID)', respName, medTerm);
            fm_nogroup = sprintf('%s ~ Condition + Age_c + %s + (1 + Condition | SubjectID)', respName, medTerm);
            fm_nocond = sprintf('%s ~ Group + Age_c + %s + (1 + Condition | SubjectID)', respName, medTerm);

            m_full_int = fitlme(T, fm_full_int, opts{:});
            m_full_noint  = fitlme(T, fm_full_noint, opts{:});
            m_nogroup = fitlme(T, fm_nogroup, opts{:});
            m_nocond = fitlme(T, fm_nocond,  opts{:});

        catch
            fm_full_int = sprintf('%s ~ Group*Condition + Age_c + %s + (1 | SubjectID)', respName, medTerm);
            fm_full_noint = sprintf('%s ~ Group + Condition + Age_c + %s + (1 | SubjectID)', respName, medTerm);
            fm_nogroup = sprintf('%s ~ Condition + Age_c + %s + (1 | SubjectID)', respName, medTerm);
            fm_nocond = sprintf('%s ~ Group + Age_c + %s + (1 | SubjectID)', respName, medTerm);

            m_full_int = fitlme(T, fm_full_int, opts{:});
            m_full_noint = fitlme(T, fm_full_noint, opts{:});
            m_nogroup = fitlme(T, fm_nogroup, opts{:});
            m_nocond = fitlme(T, fm_nocond, opts{:});
        end

        % Omnibus LRTs

        test_int = compare(m_full_noint, m_full_int);
        test_grp = compare(m_nogroup, m_full_noint);
        test_cond = compare(m_nocond, m_full_noint);

        p_int(row) = test_int.pValue(2);
        p_grp(row) = test_grp.pValue(2);
        p_cond(row) = test_cond.pValue(2);

        % Store base info

        R = struct();
        R.dataTable = T;
        R.transform = struct('respName', respName, 'epsC', epsC, 'cOcc', cOcc, 'winSec', wS, 'winSamp', wN);
        R.models.full_with_interaction = m_full_int;
        R.models.no_interaction = m_full_noint;
        R.tests.interaction = table(test_int.pValue(2),  'VariableNames', {'p'});
        R.tests.group = table(test_grp.pValue(2),  'VariableNames', {'p'});
        R.tests.condition = table(test_cond.pValue(2), 'VariableNames', {'p'});

        OUT.(measure){e, m} = R;
    end
end

% -------- BH–FDR across the 21 ERP×Microstate models (this measure) --------

p_int_fdr = bh_adjust(p_int);
p_grp_fdr = bh_adjust(p_grp);
p_cond_fdr = bh_adjust(p_cond);

% -------- EMMs and Pairwise contrasts (no margins dependency) --------

for r = 1: 21
    e = map(r, 1);
    m = map(r, 2);
    ERPname = ERPs{e};
    R = OUT.(measure){e,m};
    T = R.dataTable;

    R.tests.interaction.pFDR = p_int_fdr(r);
    R.tests.group.pFDR = p_grp_fdr(r);
    R.tests.condition.pFDR = p_cond_fdr(r);

    % Predict EMMs at Age_c = 0 and meds = 0 (no medication)

    atSpec = struct('Age_c', 0);

    % Discover which med covariates were present in this panel (robust to dropped cols)

    medNames_here = intersect({'AD', 'AP', 'MS', 'ANX', 'OTHER', 'poly_count'}, T.Properties.VariableNames, 'stable');

    for k = 1: numel(medNames_here)
        atSpec.(medNames_here{k}) = 0;
    end

    % Build med term string for reference-coded models as well

    medTerm = strjoin(medNames_here, ' + ');

    if isempty(medTerm)
        medTerm = '1 - 1';
    end % harmless placeholder if none present

    if R.tests.interaction.pFDR < alphaFDR

        % ----- Interaction present: simple-effects EMMs + pairwise -----
        R.model_used = 'full_with_interaction';

        % EMMs by strata
        R.emm.Group_by_Condition = emm_lme(R.models.full_with_interaction, {'Group','Condition'}, T, atSpec, 0.05);
        R.emm.Condition_by_Group = emm_lme(R.models.full_with_interaction, {'Condition','Group'}, T, atSpec, 0.05);
        R.emm.Group_by_Condition = addBT_emm(R.emm.Group_by_Condition, measure, R.transform);
        R.emm.Condition_by_Group = addBT_emm(R.emm.Condition_by_Group, measure, R.transform);

        % Reference-coded model for contrasts (same fixed effects + meds)

        m_full_int_ref = fitlme(T, sprintf('%s ~ Group*Condition + Age_c + %s + (1 | SubjectID)', R.transform.respName, medTerm), 'DummyVarCoding','reference','FitMethod','ML');

        % Simple-effects pairwise

        [PG, famG] = pairwise_simple(m_full_int_ref, T, measure, 'Group', 'Condition'); % Groups within each Condition
        PG.pFDR = within_families_fdr(PG.pValue, famG);
        R.pairwise.Group_within_Condition = PG;

        [PC, famC] = pairwise_simple(m_full_int_ref, T, measure, 'Condition', 'Group'); % Conditions within each Group
        PC.pFDR = within_families_fdr(PC.pValue, famC);
        R.pairwise.Condition_within_Group = PC;

    else

        % ----- No interaction: main-effects EMMs + pairwise if FDR-sig -----

        R.model_used = 'no_interaction';

        m_main_ref = fitlme(T, sprintf('%s ~ Group + Condition + Age_c + %s + (1 | SubjectID)', R.transform.respName, medTerm), 'DummyVarCoding','reference','FitMethod','REML');

        if R.tests.group.pFDR < alphaFDR
            R.emm.Group = emm_lme(m_main_ref, {'Group'}, T, atSpec, 0.05);
            R.emm.Group = addBT_emm(R.emm.Group, measure, R.transform);

            PG = pairwise_main(m_main_ref, T, measure, 'Group');
            PG.pFDR = bh_adjust(PG.pValue);
            R.pairwise.Group = PG;
        end

        if R.tests.condition.pFDR < alphaFDR
            R.emm.Condition = emm_lme(m_main_ref, {'Condition'}, T, atSpec, 0.05);
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
OUT.info.alphaFDR = alphaFDR;
OUT.info.medspec = medSpec;
OUT.info.notes = 'EMMs via predict(Conditional=false); pairwise by linear contrasts on reference-coded LMEs; back-transformed reporting; medication-adjusted.';

end % ================= END MAIN =================

% ======================= LOCAL HELPERS =======================

function [MedTbl, medNames] = meds_to_table(MedsCell, medSpec)

% Convert a cell array of medication structs to a MATLAB table of covariates.
% medSpec = 'classes' -> AD, AP, MS, ANX, OTHER (0/1 each)
% medSpec = 'poly' -> poly_count (0..5) only

n = numel(MedsCell);

% Safe getters that tolerate empty cells / missing fields

getf = @(i, f, def) (~isempty(MedsCell{i}) && isfield(MedsCell{i}, f) ) .* double(MedsCell{i}.(f)) + (isempty(MedsCell{i}) || ~isfield(MedsCell{i}, f) ) .* double(def);

switch medSpec
    case 'classes'

        AD = zeros(n, 1);
        AP = zeros(n, 1);
        MS  = zeros(n, 1);
        ANX = zeros(n, 1);
        OTH = zeros(n, 1);

        for i = 1: n
            AD(i) = getf(i, 'AD', 0);
            AP(i) = getf(i, 'AP', 0);
            MS(i) = getf(i, 'MS', 0);
            ANX(i)= getf(i, 'ANX', 0);
            OTH(i) = getf(i, 'OTHER', 0);
        end

        MedTbl = table(AD, AP, MS, ANX, OTH, 'VariableNames', {'AD', 'AP', 'MS', 'ANX', 'OTHER'});
        medNames = {'AD', 'AP', 'MS', 'ANX', 'OTHER'};

    case 'poly'
        PC = zeros(n, 1);

        for i = 1: n
            PC(i) = getf(i, 'poly_count', 0);
        end

        MedTbl = table(PC, 'VariableNames', {'poly_count'});
        medNames = {'poly_count'};
end
end

function tbl = emm_lme(M, factors, Tfit, at, alpha)

% Build a grid over "factors" and predict fixed-effects means (CIs).

if nargin < 5 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 4 || isempty(at)
    at = struct;
end

% Factor grid

nF = numel(factors);
levels = cell(nF, 1);
sz = zeros(1, nF);

for i = 1:nF
    levels{i} = categories(Tfit.(factors{i}));
    sz(i) = numel(levels{i});
end

subs = cell(1, nF);
[subs{:}] = ndgrid(levels{:});
N = numel(subs{1});

newT = table;

for i = 1 nF
    lv = string(subs{i}(:));
    newT.(factors{i}) = categorical(lv, categories(Tfit.(factors{i})));
end

newT.SubjectID = repmat(Tfit.SubjectID(1), N, 1); % any valid ID

% include other predictors at specified values

vars = setdiff(M.PredictorNames, [factors, {'SubjectID'}]);

for i = 1:numel(vars)
    vn = vars{i};
    if isfield(at, vn)
        val = at.(vn);
    else
        if isnumeric(Tfit.(vn))
            val = 0;
        else
            val = Tfit.(vn)(1);
        end
    end
    newT.(vn) = repmat(val, N, 1);
end

[yhat, yCI] = predict(M, newT, 'Conditional', false, 'Prediction', 'curve', 'Alpha', alpha);

tbl = newT(:, factors);
tbl.Estimate = yhat;
tbl.Lower = yCI(:, 1);
tbl.Upper = yCI(:, 2);
end

function T = addBT_emm(T, measure, tr)

switch measure
    case 'Duration'
        T.EMM_BT = T.Estimate;
        T.EMM_BT_Low = T.Lower;
        T.EMM_BT_High = T.Upper;
        T.Units = repmat("ms", height(T), 1);

    case 'Coverage'
        invlogit = @(x) 1./(1 + exp(-x));
        T.EMM_BT = 100*invlogit(T.Estimate);
        T.EMM_BT_Low = 100*invlogit(T.Lower);
        T.EMM_BT_High = 100*invlogit(T.Upper);
        T.Units = repmat("%", height(T), 1);

    case 'Occurrence'
        bt = @(x) max(0, exp(x) - tr.cOcc);
        T.EMM_BT = bt(T.Estimate);
        T.EMM_BT_Low = bt(T.Lower);
        T.EMM_BT_High = bt(T.Upper);
        T.Units = repmat("Hz", height(T), 1);
end
end

function [P, famIdx] = pairwise_main(Mref, T, measure, factorName)

% All pairs for a single factor in a no-interaction model (reference coding)

lv = categories(T.(factorName));
pairs = nchoosek(1:numel(lv), 2);
k = size(pairs, 1);

beta = Mref.Coefficients.Estimate;
Names = Mref.Coefficients.Name;

P = table(strings(k, 1), strings(k, 1), nan(k, 1), nan(k, 1), nan(k, 1), nan(k, 1), 'VariableNames', {'Level1', 'Level2', 'Estimate', 'SE', 'DF', 'pValue'});

for i = 1: k
    a = lv{pairs(i, 1)};
    b = lv{pairs(i, 2)};
    L = contrast_row_main(Names, factorName, a, b); % 1×p
    est = L * beta;
    [p, F, ~, df2] = coefTest(Mref, L, 0);
    se = abs(est) / max(sqrt(F), eps);
    P.Level1(i) = string(a);
    P.Level2(i) = string(b);
    P.Estimate(i) = est;
    P.SE(i) = se;
    P.DF(i) = df2;
    P.pValue(i) = p;
end

P = finalize_pairwise_table(P, measure);
famIdx = ones(height(P), 1);
end

function [P, famIdx] = pairwise_simple(Mref, T, measure, withinF, byF)

% Simple-effects pairwise within each level of byF (reference coding)

lvW = categories(T.(withinF));
lvB = categories(T.(byF));

beta = Mref.Coefficients.Estimate;
Names = Mref.Coefficients.Name;

rows = {};

for jb = 1: numel(lvB)
    pr = nchoosek(1:numel(lvW), 2);

    for i = 1: size(pr, 1)
        rows(end + 1, :) = {lvB{jb}, lvW{pr(i, 1)}, lvW{pr(i, 2)}};
    end
end

n = size(rows, 1);
P = table(strings(n, 1), strings(n, 1), strings(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), 'VariableNames', {byF, 'Level1', 'Level2', 'Estimate', 'SE', 'DF', 'pValue'});

for i = 1: n
    bf = rows{i, 1};
    a = rows{i, 2};
    b = rows{i, 3};
    L = contrast_row_simple(Names, withinF, byF, a, b, bf); % 1×p
    est = L * beta;
    [p, F, ~, df2] = coefTest(Mref, L, 0);
    se = abs(est) / max(sqrt(F), eps);
    P.(byF)(i) = string(bf);
    P.Level1(i) = string(a);
    P.Level2(i) = string(b);
    P.Estimate(i)= est;
    P.SE(i) = se;
    P.DF(i) = df2;
    P.pValue(i) = p;
end

P = finalize_pairwise_table(P, measure);

% family index per each byF level

block = nchoosek(1: numel(lvW), 2);
per = size(block, 1);
famIdx = zeros(height(P), 1);

for jb = 1: numel(lvB)
    idx = (1: per) + (jb - 1) * per;
    famIdx(idx) = jb;
end
end

function P = finalize_pairwise_table(P, measure)

% Add measure-specific effect-size columns + CIs from est±t*SE

crit = tinv(0.975, P.DF);
lo = P.Estimate - crit.*P.SE;
hi = P.Estimate + crit.*P.SE;

switch measure
    case 'Duration'
        P.Delta = P.Estimate;
        P.Delta_Low = lo;
        P.Delta_High = hi;

    case 'Coverage'
        P.OR = exp(P.Estimate);
        P.OR_Low = exp(lo);
        P.OR_High = exp(hi);

    case 'Occurrence'
        P.RR = exp(P.Estimate);
        P.RR_Low = exp(lo);
        P.RR_High = exp(hi);
end
end

function L = contrast_row_main(Names, factor, L1, L2)

% No-interaction model, reference coding.

p = numel(Names);
L = zeros(1,p);

    function addcoef(nm, sgn)
        j = strcmp(Names, nm);

        if any(j)
            L(j) = L(j) + sgn;
        end
    end

addcoef(sprintf('%s_%s', factor, L1), + 1);
addcoef(sprintf('%s_%s', factor, L2), - 1);
end

function L = contrast_row_simple(Names, withinF, byF, L1, L2, BY)

% Interaction model, reference coding. Simple effect (L1-L2) of withinF at BY.

p = numel(Names);
L = zeros(1, p);

    function addcoef(nm, sgn)

        j = strcmp(Names, nm);

        if any(j)
            L(j) = L(j) + sgn;
        end
    end

    function addint(f1, lv1, f2, lv2, sgn)
        nm1 = sprintf('%s_%s:%s_%s', f1, lv1, f2, lv2);
        nm2 = sprintf('%s_%s:%s_%s', f2, lv2, f1, lv1);
        j = strcmp(Names, nm1) | strcmp(Names, nm2);

        if any(j)
            L(j) = L(j) + sgn;
        end
    end

addcoef(sprintf('%s_%s', withinF, L1), +1);
addcoef(sprintf('%s_%s', withinF, L2), -1);
addint(withinF, L1, byF, BY, +1);
addint(withinF, L2, byF, BY, -1);
end

function p_adj = within_families_fdr(p, famIdx)

p_adj = nan(size(p));

for u = unique(famIdx(:))'
    ii = famIdx == u;

    if any(ii)
        p_adj(ii) = bh_adjust(p(ii));
    end

end
end

function p_adj = bh_adjust(p)
p = p(:);
m = numel(p);
[ps, idx] = sort(p);
q = ps .* m ./ (1:m)';

for i = m - 1: -1: 1
    q(i) = min(q(i), q(i+1));
end

p_adj = nan(size(p));
p_adj(idx) = q;
end