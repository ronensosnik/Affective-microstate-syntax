function print_summary_measure(OUT_measure)

% PRINT_SUMMARY_MEASURE
% Pretty on-screen summary per ERP × Microstate for a single measure.
% Shows EMMs (back-transformed) and pairwise contrasts with directionality.
%
% Usage: print_summary_measure(Out_Coverage)

ERPs = OUT_measure.info.ERPs;
Micros = OUT_measure.info.Microstates;
measure = OUT_measure.info.measure;

if isfield(OUT_measure.info, 'alphaFDR')
    alphaFDR = OUT_measure.info.alphaFDR;
else
    alphaFDR = 0.05;
end

barline = repmat('-', 1, 78);
fprintf('\n%s\nSUMMARY (%s) — ERP BH-FDR alpha = %.3f\n%s\n', barline, measure, alphaFDR, barline);

for e = 1: numel(ERPs)
    for m = 1: numel(Micros)

        R = OUT_measure.(measure){e, m};

        if isempty(R)
            continue;
        end

        p_interaction = get_nested(R, {'tests', 'interaction', 'pFDR'});
        p_Group = get_nested(R, {'tests', 'group', 'pFDR'});
        p_Condition = get_nested(R, {'tests', 'condition', 'pFDR'});

        [M, modelLabel] = pick_model(R);

        header = sprintf('%s | %s × %s', measure, ERPs{e}, Micros{m});

        fprintf('\n%s\n%s\n', header, repmat('=', 1, numel(header)));
        fprintf('Model used: %s\n', modelLabel);

        fprintf('Omnibus (ERP-wise BH-FDR): Interaction p = %s %s | Group p = %s %s | Condition p = %s %s\n', fp(p_interaction), star(isfinite(p_interaction) && p_interaction < alphaFDR), fp(p_Group), star(isfinite(p_Group) && p_Group < alphaFDR), fp(p_Condition), star(isfinite(p_Condition) && p_Condition < alphaFDR));

        ageLine = age_effect_line(M);
        if strlength(ageLine) > 0
            fprintf('%s\n', ageLine);
        end

        if strcmpi(measure, 'Num_occurrence') && isfield(R, 'transform') && isfield(R.transform, 'cOcc')
            fprintf('Note: Occurrence modeled as log(rate + c), where c = 0.5 / window-seconds = %.4f Hz.\n', R.transform.cOcc);
        end

        if isfinite(p_interaction) && p_interaction < alphaFDR

            if isfield(R, 'emm') && isfield(R.emm, 'Group_by_Condition')
                T = R.emm.Group_by_Condition;
                fprintf('\nEMMs: Group within each Condition (%s)\n', units_for(measure));
                disp_cols(T, {'Group', 'Condition', 'Estimate_BT', 'CI_lo_BT', 'CI_hi_BT', 'Units'});
            end

            if isfield(R, 'emm') && isfield(R.emm, 'Condition_by_Group')
                T = R.emm.Condition_by_Group;
                fprintf('\nEMMs: Condition within each Group (%s)\n', units_for(measure));
                disp_cols(T, {'Condition', 'Group', 'Estimate_BT', 'CI_lo_BT', 'CI_hi_BT', 'Units'});
            end

            if isfield(R, 'pairwise') && isfield(R.pairwise, 'Group_within_Condition') && isfield(R, 'emm') && isfield(R.emm, 'Group_by_Condition')
                fprintf('\nPairwise: Groups within each Condition (direction and CI)\n');
                print_pairs_simple(R.pairwise.Group_within_Condition, R.emm.Group_by_Condition, measure, 'Group', 'Condition', alphaFDR);
            end

            if isfield(R, 'pairwise') && isfield(R.pairwise, 'Condition_within_Group') && isfield(R, 'emm') && isfield(R.emm, 'Condition_by_Group')
                fprintf('\nPairwise: Conditions within each Group (direction and CI)\n');
                print_pairs_simple(R.pairwise.Condition_within_Group, R.emm.Condition_by_Group, measure, 'Condition', 'Group', alphaFDR);
            end

        else

            if isfinite(p_Group) && p_Group < alphaFDR && isfield(R, 'emm') && isfield(R.emm, 'Group')
                fprintf('\nEMMs: Group (%s)\n', units_for(measure));
                disp_cols(R.emm.Group, {'Group', 'Estimate_BT', 'CI_lo_BT', 'CI_hi_BT', 'Units'});
            end

            if isfinite(p_Condition) && p_Condition < alphaFDR && isfield(R, 'emm') && isfield(R.emm, 'Condition')
                fprintf('\nEMMs: Condition (%s)\n', units_for(measure));
                disp_cols(R.emm.Condition, {'Condition', 'Estimate_BT', 'CI_lo_BT', 'CI_hi_BT', 'Units'});
            end

            if isfinite(p_Group) && p_Group < alphaFDR && isfield(R, 'pairwise') && isfield(R.pairwise, 'Group') && isfield(R, 'emm') && isfield(R.emm, 'Group')
                fprintf('\nPairwise: Group (direction and CI)\n');
                print_pairs_main(R.pairwise.Group, R.emm.Group, measure, 'Group', alphaFDR);
            end

            if isfinite(p_Condition) && p_Condition < alphaFDR && isfield(R, 'pairwise') && isfield(R.pairwise, 'Condition') && isfield(R, 'emm') && isfield(R.emm, 'Condition')
                fprintf('\nPairwise: Condition (direction and CI)\n');
                print_pairs_main(R.pairwise.Condition, R.emm.Condition, measure, 'Condition', alphaFDR);
            end

        end

    end
end

fprintf('%s\n', barline);

end


function [M, label] = pick_model(R)

M = [];
label = "";

if isfield(R, 'model_used') && strcmpi(R.model_used, 'full_with_interaction') && isfield(R, 'models') && isfield(R.models, 'full_with_interaction')
    M = R.models.full_with_interaction;
    label = "full with interaction";
    return;
end

if isfield(R, 'models') && isfield(R.models, 'no_interaction_REML')
    M = R.models.no_interaction_REML;
    label = "no interaction (REML)";
    return;
end

if isfield(R, 'models') && isfield(R.models, 'no_interaction_ML')
    M = R.models.no_interaction_ML;
    label = "no interaction (ML)";
    return;
end

if isfield(R, 'models') && isfield(R.models, 'no_interaction')
    M = R.models.no_interaction;
    label = "no interaction";
    return;
end

label = "unknown";

end


function s = age_effect_line(M)

s = "";

if isempty(M)
    return;
end

try
    C = M.Coefficients;
    idx = strcmp(C.Name, 'Age_c');
    if any(idx)
        b = C.Estimate(idx);
        p = C.pValue(idx);
        s = sprintf('Covariate: Age_c beta = %s, p = %s (uncorrected).', fnum(b, 4), fp(p));
    end
catch
    s = "";
end

end


function disp_cols(T, cols)

have = ismember(cols, T.Properties.VariableNames);
cols = cols(have);

if isempty(cols)
    disp(T);
else
    disp(T(:, cols));
end

end


function print_pairs_main(P, T, measure, factorName, alphaFDR)

[pAdj, pAdjName] = pick_padj(P);

mask = isfinite(pAdj) & (pAdj < alphaFDR);

if ~any(mask)
    fprintf('No pairwise survived FDR.\n');
    return;
end

u = units_for(measure);
rows = find(mask);

for r = rows(:).'

    [L1, L2] = pick_levels(P, r, factorName);
    E1 = pick_emm(T, factorName, L1);
    E2 = pick_emm(T, factorName, L2);

    eff = effect_text(P, r, measure);
    dirTxt = direction_text(E1, E2, L1, L2, u);

    pRaw = pick_praw(P, r);

    fprintf('  %s vs %s: %s; EMMs: %s %s (%s) vs %s %s (%s); p = %s (%s = %s)\n', L1, L2, eff, fnum(E1, 2), u, L1, fnum(E2, 2), u, L2, fp(pRaw), pAdjName, fp(pAdj(r)));

    if strlength(dirTxt) > 0
        fprintf('    Direction: %s\n', dirTxt);
    end

    fprintf('\n');
end

end


function print_pairs_simple(P, T, measure, withinF, byF, alphaFDR)

[pAdj, pAdjName] = pick_padj(P);
sig = isfinite(pAdj) & pAdj < alphaFDR;

if ~any(sig)
    fprintf('  No pairwise survived FDR.\n');
    return;
end

u = units_for(measure);
rows = find(sig);

for r = rows(:).'

    ctx = string(P.(byF)(r));
    [L1, L2] = pick_levels(P, r, withinF);

    if strcmpi(withinF, 'Group')
        E1 = pick_emm2(T, 'Group', L1, 'Condition', ctx);
        E2 = pick_emm2(T, 'Group', L2, 'Condition', ctx);
    else
        E1 = pick_emm2(T, 'Condition', L1, 'Group', ctx);
        E2 = pick_emm2(T, 'Condition', L2, 'Group', ctx);
    end

    eff = effect_text(P, r, measure);
    dirTxt = direction_text(E1, E2, L1, L2, u);
    pRaw = pick_praw(P, r);

    fprintf('%s vs %s @ %s = %s: %s; EMMs: %s %s (%s) vs %s %s (%s); p = %s (%s = %s)\n', L1, L2, byF, ctx, eff, fnum(E1, 2), u, L1, fnum(E2, 2), u, L2, fp(pRaw), pAdjName, fp(pAdj(r)));

    if strlength(dirTxt) > 0
        fprintf('    Direction: %s\n', dirTxt);
    end

end

end


function [pAdj, pAdjName] = pick_padj(P)

if any(strcmp(P.Properties.VariableNames, 'pFDR'))
    pAdj = P.pFDR;
    pAdjName = "pFDR";
elseif any(strcmp(P.Properties.VariableNames, 'pAdj'))
    pAdj = P.pAdj;
    pAdjName = "pAdj";
elseif any(strcmp(P.Properties.VariableNames, 'qValue'))
    pAdj = P.qValue;
    pAdjName = "qValue";
else
    pAdj = P.pValue;
    pAdjName = "p";
end

end


function x = pick_emm(T, fac, lvl)

mask = T.(fac) == categorical(lvl, categories(T.(fac)));
x = NaN;

if any(mask)
    col = pick_emm_col(T);
    x = T.(col)(find(mask, 1, 'first'));
end

end


function x = pick_emm2(T, f1, l1, f2, l2)

mask = T.(f1) == categorical(l1, categories(T.(f1))) & T.(f2) == categorical(l2, categories(T.(f2)));
x = NaN;

if any(mask)
    col = pick_emm_col(T);
    x = T.(col)(find(mask, 1, 'first'));
end

end


function col = pick_emm_col(T)

if any(strcmp(T.Properties.VariableNames, 'EMM_BT'))
    col = 'EMM_BT';
elseif any(strcmp(T.Properties.VariableNames, 'Estimate_BT'))
    col = 'Estimate_BT';
elseif any(strcmp(T.Properties.VariableNames, 'Estimate'))
    col = 'Estimate';
else
    col = T.Properties.VariableNames{1};
end

end


function p = pick_praw(P, r)

p = P.pValue(r);

end


function [L1, L2] = pick_levels(P, r, contrastVar)

L1 = "";
L2 = "";

if all(ismember({'Level1', 'Level2'}, P.Properties.VariableNames))
    L1 = string(P.Level1(r));
    L2 = string(P.Level2(r));
    return;
end

if any(strcmp(P.Properties.VariableNames, contrastVar))
    s = string(P.(contrastVar)(r));
elseif any(strcmp(P.Properties.VariableNames, 'Contrast'))
    s = string(P.Contrast(r));
else
    s = "";
    for j = 1: width(P)
        if iscellstr(P{:, j}) || isstring(P{:, j}) || iscategorical(P{:, j})
            s = string(P{r, j});
            break;
        end
    end
end

[L1, L2] = split_contrast(s);

end


function [a, b] = split_contrast(s)

s = strtrim(string(s));

s = replace(s, char(8211), '-');
s = replace(s, char(8212), '-');

% Require at least one space around the separator so labels like "BD-I" are not split.
parts = regexp(s, '\s+(?:-|vs\.?|VS\.?|Vs\.?|→)\s+', 'split');
parts = parts(~cellfun(@isempty, parts));

if numel(parts) >= 2
    a = strtrim(string(parts{1}));
    b = strtrim(string(parts{2}));
else
    a = s;
    b = "";
end

end


function eff = effect_text(P, r, measure)

meas = lower(string(measure));

if any(strcmp(P.Properties.VariableNames, 'Estimate'))
    est = P.Estimate(r);
elseif any(strcmp(P.Properties.VariableNames, 'Estimate_BT'))
    est = P.Estimate_BT(r);
else
    est = NaN;
end

lo_est = NaN;
hi_est = NaN;

hasLo = any(strcmp(P.Properties.VariableNames, 'Lower'));
hasHi = any(strcmp(P.Properties.VariableNames, 'Upper'));

if hasLo
    lo_est = P.Lower(r);
end
if hasHi
    hi_est = P.Upper(r);
end

if any(strcmp(P.Properties.VariableNames, 'SE'))
    se = P.SE(r);
else
    se = NaN;
end

if any(strcmp(P.Properties.VariableNames, 'DF'))
    df = P.DF(r);
else
    df = NaN;
end

if ~isfinite(lo_est) || ~isfinite(hi_est)
    if isfinite(se)
        if isfinite(df)
            crit = tinv(0.975, df);
        else
            crit = 1.96;
        end
        lo_est = est - crit .* se;
        hi_est = est + crit .* se;
    end
end

if meas == "duration"
    d = est;
    lo = lo_est;
    hi = hi_est;
    eff = sprintf('Δ = %s ms [%s, %s]', sgn(d), fnum(lo, 2), fnum(hi, 2));
elseif meas == "coverage"
    rr = exp(est);
    lo = exp(lo_est);
    hi = exp(hi_est);
    eff = sprintf('OR = %s [%s, %s]', fnum(rr, 2), fnum(lo, 2), fnum(hi, 2));
else
    rr = exp(est);
    lo = exp(lo_est);
    hi = exp(hi_est);
    eff = sprintf('RR = %s [%s, %s]', fnum(rr, 2), fnum(lo, 2), fnum(hi, 2));
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


function u = units_for(measure)

switch measure
    case 'Duration'
        u = 'ms';
    case 'Coverage'
        u = '%';
    otherwise
        u = 'Hz';
end

end


function s = direction_text(x1, x2, l1, l2, u)

if isnan(x1) || isnan(x2)
    s = "";
    return;
end

if x1 > x2
    s = sprintf('%s > %s by %s %s', l1, l2, fnum(x1 - x2, 2), u);
elseif x2 > x1
    s = sprintf('%s > %s by %s %s', l2, l1, fnum(x2 - x1, 2), u);
else
    s = sprintf('%s = %s', l1, l2);
end

end


function s = fp(p)

if isnan(p)
    s = 'NA';
elseif p < 1e-4
    s = '<1e-4';
else
    s = sprintf('%.3f', p);
end

end


function s = fnum(x, d)

if isnan(x)
    s = 'NA';
else
    s = sprintf(['%0.', num2str(d), 'f'], x);
end

end


function s = sgn(x)

if isnan(x)
    s = 'NA';
elseif x >= 0
    s = sprintf('+%.2f', x);
else
    s = sprintf('%.2f', x);
end

end


function s = star(tf)

if tf
    s = '★';
else
    s = '';
end

end


function p_adj = bh_adjust(p)

p = p(:);
keep = isfinite(p);

p_adj = nan(size(p));

if ~any(keep)
    return;
end

pp = p(keep);
[ps, idx] = sort(pp, 'ascend');

m = numel(ps);
q = ps .* m ./ (1: m)';

for i = m - 1: -1: 1
    q(i) = min(q(i), q(i + 1));
end

q(q > 1) = 1;

tmp = nan(size(pp));
tmp(idx) = q;

p_adj(keep) = tmp;

end
