function [overall, byERP] = summarize_q1_robustness(Base, Sens, ERPs, Groups, Conds)

% SUMMARIZE_Q1_ROBUSTNESS
% Compare sensitivity pipeline vs baseline for Q1 (within-cell deviation from independence).
%
% Primary:
%   - Among BASELINE-significant edges: direction retention and significance retention.
%
% Added:
%   - Jaccard overlap of significant edge sets (per cell; summarized overall and by ERP)
%   - Lost/gained significant edge counts (per cell; summarized)
%   - Percent agreement and Cohen's kappa for edge significance (per cell; summarized)
%   - Effect-size concordance (Spearman rho of log(IRR), RMSE, MAD; per cell; summarized)
%   - Separate summaries for:
%       (a) All cells (includes empty–empty panels; Jaccard(empty,empty)=1)
%       (b) Non-empty cells only (union of sig edges > 0)

alpha = 0.05; %#ok<NASGU>  % retained for clarity / future extensions

edgeMask = ~eye(7);

total_base = 0;
retain_dir = 0;
retain_sig = 0;

cellStats = []; % accumulate per ERP×Group×Cond cell stats
per_erp = struct();

for e = 1: numel(ERPs)

    erp = ERPs{e};

    bE = 0;
    dE = 0;
    sE = 0;

    cellRows = [];

    for g = 1: numel(Groups)

        grp = Groups{g};

        for c = 1: numel(Conds)

            cond = Conds{c};

            B = safeget(Base, [], erp, grp, cond);
            S = safeget(Sens, [], erp, grp, cond);

            if isempty(B) || isempty(S)
                continue;
            end

            if ~isfield(B, 'sig') || ~isfield(S, 'sig')
                continue;
            end

            if ~isfield(B, 'IRR') || ~isfield(S, 'IRR')
                continue;
            end

            sigB = logical(B.sig) & edgeMask;
            sigS = logical(S.sig) & edgeMask;

            % ---- Baseline significant edges: retain direction + retain significance ----

            if any(sigB(:))

                idx = find(sigB);

                for k = 1: numel(idx)

                    total_base = total_base + 1;
                    bE = bE + 1;

                    [i, j] = ind2sub([7 7], idx(k));

                    dirB = dir_from_ratio(B.IRR(i, j));
                    dirS = dir_from_ratio(S.IRR(i, j));

                    if ~isnan(dirB) && ~isnan(dirS) && dirB == dirS
                        retain_dir = retain_dir + 1;
                        dE = dE + 1;
                    end

                    if sigS(i, j)
                        retain_sig = retain_sig + 1;
                        sE = sE + 1;
                    end

                end
            end

            % ---- Added per-cell robustness stats ----

            [jacc, lostN, gainedN, agreePct, kappaVal, unionN] = sig_overlap_stats(sigB, sigS);
            [rho, rmseLog, madLog] = effect_concordance(B.IRR, S.IRR, edgeMask);

            cellRows = [cellRows; struct( ...
                'ERP', string(erp), ...
                'Group', string(grp), ...
                'Condition', string(cond), ...
                'Jaccard', jacc, ...
                'Kappa', kappaVal, ...
                'AgreePct', agreePct, ...
                'LostEdges', lostN, ...
                'GainedEdges', gainedN, ...
                'UnionEdges', unionN, ...
                'NonEmpty', unionN > 0, ...
                'RhoLog', rho, ...
                'RMSELog', rmseLog, ...
                'MADLog', madLog ...
                )];

        end
    end

    if isempty(cellRows)
        per_erp.(erp) = struct( ...
            'BaselineEdges', bE, ...
            'DirectionRetained', dE, ...
            'SignificanceRetained', sE, ...
            'DirPct', pct(dE, bE), ...
            'SigPct', pct(sE, bE), ...
            'CellsCompared', 0, ...
            'Jaccard_Median_AllCells', NaN, ...
            'Kappa_Median_AllCells', NaN, ...
            'NonEmptyCells', 0, ...
            'Jaccard_Median_NonEmpty', NaN, ...
            'Kappa_Median_NonEmpty', NaN, ...
            'KappaDefined_NonEmptyCells', 0, ...
            'RhoLog_Median', NaN, ...
            'RMSELog_Median', NaN, ...
            'MADLog_Median', NaN, ...
            'LostEdges_Total', 0, ...
            'GainedEdges_Total', 0 ...
            );
    else
        Tcell = struct2table(cellRows);
        cellStats = [cellStats; Tcell];

        nonEmptyMask = logical(Tcell.NonEmpty);
        kappaDefinedNonEmpty = sum(nonEmptyMask & isfinite(Tcell.Kappa));

        per_erp.(erp) = struct( ...
            'BaselineEdges', bE, ...
            'DirectionRetained', dE, ...
            'SignificanceRetained', sE, ...
            'DirPct', pct(dE, bE), ...
            'SigPct', pct(sE, bE), ...
            'CellsCompared', height(Tcell), ...
            'Jaccard_Median_AllCells', median(Tcell.Jaccard, 'omitnan'), ...
            'Kappa_Median_AllCells', median(Tcell.Kappa, 'omitnan'), ...
            'NonEmptyCells', sum(nonEmptyMask), ...
            'Jaccard_Median_NonEmpty', median(Tcell.Jaccard(nonEmptyMask), 'omitnan'), ...
            'Kappa_Median_NonEmpty', median(Tcell.Kappa(nonEmptyMask), 'omitnan'), ...
            'KappaDefined_NonEmptyCells', kappaDefinedNonEmpty, ...
            'RhoLog_Median', median(Tcell.RhoLog, 'omitnan'), ...
            'RMSELog_Median', median(Tcell.RMSELog, 'omitnan'), ...
            'MADLog_Median', median(Tcell.MADLog, 'omitnan'), ...
            'LostEdges_Total', sum(Tcell.LostEdges, 'omitnan'), ...
            'GainedEdges_Total', sum(Tcell.GainedEdges, 'omitnan') ...
            );
    end

end

% ---- Overall summary ----

overall = struct();
overall.BaselineEdges = total_base;
overall.DirectionRetained = retain_dir;
overall.SignificanceRetained = retain_sig;
overall.DirPct = pct(retain_dir, total_base);
overall.SigPct = pct(retain_sig, total_base);

if ~isempty(cellStats)

    nonEmptyMask = logical(cellStats.NonEmpty);
    kappaDefinedNonEmpty = sum(nonEmptyMask & isfinite(cellStats.Kappa));

    overall.CellsCompared = height(cellStats);

    overall.Jaccard_Median_AllCells = median(cellStats.Jaccard, 'omitnan');
    overall.Kappa_Median_AllCells = median(cellStats.Kappa, 'omitnan');

    overall.NonEmptyCells = sum(nonEmptyMask);

    overall.Jaccard_Median_NonEmpty = median(cellStats.Jaccard(nonEmptyMask), 'omitnan');
    overall.Kappa_Median_NonEmpty = median(cellStats.Kappa(nonEmptyMask), 'omitnan');
    overall.KappaDefined_NonEmptyCells = kappaDefinedNonEmpty;

    overall.RhoLog_Median = median(cellStats.RhoLog, 'omitnan');
    overall.RMSELog_Median = median(cellStats.RMSELog, 'omitnan');
    overall.MADLog_Median = median(cellStats.MADLog, 'omitnan');

    overall.LostEdges_Total = sum(cellStats.LostEdges, 'omitnan');
    overall.GainedEdges_Total = sum(cellStats.GainedEdges, 'omitnan');
end

% ---- byERP table (expanded to match your Q2 format) ----

byERP = table('Size', [numel(ERPs) 14], ...
    'VariableTypes', {'string','double','double','double','string', ...
                      'double','double','double','double','double', ...
                      'double','double','double','double'}, ...
    'VariableNames', {'ERP','BaselineEdges','DirectionRetained','SignificanceRetained','Percentages', ...
                      'JaccardMed_All','KappaMed_All','JaccardMed_NonEmpty','KappaMed_NonEmpty','NonEmptyCells', ...
                      'RhoLogMed','RMSELogMed','MADLogMed','CellsCompared'});

for e = 1: numel(ERPs)

    erp = ERPs{e};
    S = per_erp.(erp);

    byERP.ERP(e) = string(erp);
    byERP.BaselineEdges(e) = S.BaselineEdges;
    byERP.DirectionRetained(e) = S.DirectionRetained;
    byERP.SignificanceRetained(e) = S.SignificanceRetained;

    byERP.Percentages(e) = sprintf('Dir=%s, Sig=%s', ...
        pct(S.DirectionRetained, S.BaselineEdges), pct(S.SignificanceRetained, S.BaselineEdges));

    byERP.JaccardMed_All(e) = safe_num(S, 'Jaccard_Median_AllCells');
    byERP.KappaMed_All(e) = safe_num(S, 'Kappa_Median_AllCells');

    byERP.JaccardMed_NonEmpty(e) = safe_num(S, 'Jaccard_Median_NonEmpty');
    byERP.KappaMed_NonEmpty(e) = safe_num(S, 'Kappa_Median_NonEmpty');
    byERP.NonEmptyCells(e) = safe_num(S, 'NonEmptyCells');

    byERP.RhoLogMed(e) = safe_num(S, 'RhoLog_Median');
    byERP.RMSELogMed(e) = safe_num(S, 'RMSELog_Median');
    byERP.MADLogMed(e) = safe_num(S, 'MADLog_Median');
    byERP.CellsCompared(e) = safe_num(S, 'CellsCompared');

end

end


% ========================= HELPERS (local) =========================

function S = safeget(S0, defaultVal, varargin)

S = defaultVal;

try
    for k = 1: numel(varargin)
        S0 = S0.(varargin{k});
    end
    S = S0;
catch
end

end


function x = safe_num(S, fieldname)

if isfield(S, fieldname)
    x = S.(fieldname);
else
    x = NaN;
end

end


function d = dir_from_ratio(r)

d = NaN;

if ~isfinite(r) || r <= 0
    return;
end

x = log(r);

if x > 0
    d = +1;
elseif x < 0
    d = -1;
else
    d = 0;
end

end


function s = pct(a, b)

if ~isfinite(a) || ~isfinite(b) || b <= 0
    s = 'NA';
    return;
end

s = sprintf('%.1f%%', 100 * a / b);

end


function [jacc, lostN, gainedN, agreePct, kappaVal, unionN] = sig_overlap_stats(sigB, sigS)

sigB = logical(sigB(:));
sigS = logical(sigS(:));

intN = sum(sigB & sigS);
unionN = sum(sigB | sigS);

if unionN == 0
    jacc = 1; % convention: empty–empty panels match perfectly
else
    jacc = intN / unionN;
end

lostN = sum(sigB & ~sigS);
gainedN = sum(~sigB & sigS);

agreePct = 100 * mean(sigB == sigS);

kappaVal = cohens_kappa(sigB, sigS);

end


function kappaVal = cohens_kappa(a, b)

a = logical(a(:));
b = logical(b(:));

po = mean(a == b);

pa1 = mean(a);
pb1 = mean(b);

pe = pa1 * pb1 + (1 - pa1) * (1 - pb1);

den = (1 - pe);

if den <= 0
    kappaVal = NaN; % undefined in degenerate cases (e.g., all-zero in both)
else
    kappaVal = (po - pe) / den;
end

end


function [rho, rmseLog, madLog] = effect_concordance(IRR_B, IRR_S, edgeMask)

rho = NaN;
rmseLog = NaN;
madLog = NaN;

try
    vB = log(IRR_B(edgeMask));
    vS = log(IRR_S(edgeMask));
catch
    return;
end

good = isfinite(vB) & isfinite(vS);

if sum(good) < 3
    return;
end

vB = vB(good);
vS = vS(good);

try
    rho = corr(vB, vS, 'Type', 'Spearman');
catch
    rho = NaN;
end

rmseLog = sqrt(mean((vB - vS).^2, 'omitnan'));
madLog = median(abs(vB - vS), 'omitnan');

end
