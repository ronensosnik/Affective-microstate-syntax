function CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, varargin)

% buildAnchoredContrasts (one-stage Poisson GLM, BH-FDR per ERP)
%
% One-stage per pair (Ref + Target):
%   y_{j,XY} ~ Poisson(mu_{j,XY})
%   log(mu_{j,XY}) = log(e_{j,XY}) + alpha0_{XY} + alpha1_{XY} * Age_ref_z(j) + delta_{XY} * 1{G_j = Target}
%
% Age_ref_z is z-scored using the reference group's mean and SD within each ERP × Condition panel.
% Multiplicity: BH-FDR pooled per ERP across all TargetGroups × Conditions × 42 edges.

p = inputParser;
addParameter(p, 'RefGroup', varargin{2});
addParameter(p, 'TargetGroups', {});
addParameter(p, 'Alpha', varargin{6});
addParameter(p, 'MaxIter', varargin{8});
addParameter(p, 'MinObs', varargin{10});
addParameter(p, 'MinExp', varargin{12});
addParameter(p, 'LogEClamp', varargin{14});
parse(p, varargin{:});
opt = p.Results;

refGrp = 'HC';
TargetList = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'Siblings'};

[srcIdx, tgtIdx] = find(~eye(7));
GLM_OPTS = statset('MaxIter', opt.MaxIter, 'TolX', 1e-10, 'TolFun', 1e-10);
MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

CompStats = struct();

% Pool Q2 p-values separately for each ERP

poolQ2 = repmat(struct('p', [], 'map', []), numel(ERPs), 1);

% map rows: [c g s t] for that ERP index (implicit by poolQ2(e))
% where g indexes TargetList, not Groups.

for e = 1: numel(ERPs)

    erp = ERPs{e};

    for c = 1: numel(Conds)

        condName = Conds{c};
        condKey = lower(condName);

        RefCell = getfield_safe(Out, {erp, condKey, refGrp});

        if isfield(AgeMap, refGrp)
            agesRef = AgeMap.(refGrp)(:);
        else
            agesRef = [];
        end

        muRef = mean(agesRef, 'omitnan');
        sdRef = std(agesRef, 'omitnan');

        if ~isfinite(muRef) || ~isfinite(sdRef) || sdRef <= eps
            muRef = 0;
            sdRef = 1;
        end

        % Prebuild reference long table for this ERP × Condition

        Table_ref = table();

        if ~isempty(RefCell) && ~isempty(agesRef)

            for sj = 1: numel(RefCell)

                S = RefCell{sj};

                if badSubject(S)
                    continue;
                end

                [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
                Dur_src_tgt = cell2mat(S.transition_duration_matrix(sub2ind([7, 7], srcIdx, tgtIdx)));
                Dur_src = Dur_src_tgt(:, 1);
                Dur_trg = Dur_src_tgt(:, 2);

                if ~any(keep)
                    continue;
                end

                eKeep = evec(keep);
                logE = log(eKeep);
                logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                ageRaw = agesRef(min(sj, numel(agesRef)));
                ageS = (ageRaw - muRef) ./ sdRef;

                R = table(yvec(keep), eKeep, logE, repmat(ageS, nnz(keep), 1), srcIdx(keep), tgtIdx(keep), zeros(nnz(keep), 1), Dur_src(keep, :), Dur_trg(keep, :), 'VariableNames', {'y', 'e', 'logE', 'AgeS', 'src', 'tgt', 'Target', 'Dur_src', 'Dur_trg'});
                Table_ref = [Table_ref; R];
            end
        end

        for g = 1: numel(TargetList)

            tgtGrp = char(TargetList{g});
            TgtCell = getfield_safe(Out, {erp, condKey, tgtGrp});

            if isfield(AgeMap, tgtGrp)
                agesTgt = AgeMap.(tgtGrp)(:);
            else
                agesTgt = [];
            end

            IRR = NaN(7);
            CIlo = NaN(7);
            CIhi = NaN(7);
            pMat = NaN(7);
            pFDR_Mat = NaN(7);
            sig = false(7);

            if isempty(Table_ref) || isempty(TgtCell) || isempty(agesTgt)
                CompStats = storeCell(CompStats, erp, condName, tgtGrp, refGrp, IRR, CIlo, CIhi, pMat, pFDR_Mat, sig, emptyEdges(), NaN);
                continue;
            end

            % Prebuild target long table

            Table_target = table();

            for sj = 1: numel(TgtCell)
                S = TgtCell{sj};

                if badSubject(S)
                    continue;
                end

                [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                if ~any(keep)
                    continue;
                end

                eKeep = evec(keep);
                logE = log(eKeep);
                logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                ageRaw = agesTgt(min(sj, numel(agesTgt)));
                ageS = (ageRaw - muRef) ./ sdRef;

                R = table(yvec(keep), eKeep, logE, repmat(ageS, nnz(keep), 1), srcIdx(keep), tgtIdx(keep), ones(nnz(keep), 1), Dur_src(keep, :), Dur_trg(keep, :), 'VariableNames', {'y', 'e', 'logE', 'AgeS', 'src', 'tgt', 'Target', 'Dur_src', 'Dur_trg'});
                Table_target = [Table_target; R];
            end

            if isempty(Table_target)
                CompStats = storeCell(CompStats, erp, condName, tgtGrp, refGrp, IRR, CIlo, CIhi, pMat, pFDR_Mat, sig, emptyEdges(), NaN);
                continue;
            end

            % Fit per pair: one-stage Poisson GLM with offset logE

            for k = 1: numel(srcIdx)
                s = srcIdx(k);
                t = tgtIdx(k);

                Rref = Table_ref(Table_ref.src == s & Table_ref.tgt == t, :);
                Rtgt = Table_target(Table_target.src == s & Table_target.tgt == t, :);

                if isempty(Rref) || isempty(Rtgt)
                    continue;
                end

                R = [Rref; Rtgt];

                if sum(R.y) < opt.MinObs || sum(R.e) < opt.MinExp
                    continue;
                end

                y_ref = sum(R.y(R.Target == 0));
                y_tgt = sum(R.y(R.Target == 1));

                if y_ref == 0 || y_tgt == 0
                    continue;
                end

                try
                    mdl = fitglm(R, 'y ~ 1 + AgeS + Target', 'Distribution', 'poisson', 'Link', 'log', 'Offset', R.logE, 'Options', GLM_OPTS);

                    idx = strcmp(mdl.CoefficientNames, 'Target');

                    if ~any(idx)
                        continue;
                    end

                    idx1 = find(idx, 1, 'first');

                    b = mdl.Coefficients.Estimate(idx1);
                    pEdge = mdl.Coefficients.pValue(idx1);
                    CI = coefCI(mdl, 0.05);

                    IRR(s, t) = exp(b);
                    CIlo(s, t) = exp(CI(idx1, 1));
                    CIhi(s, t) = exp(CI(idx1, 2));
                    pMat(s, t) = pEdge;

                    if isfinite(pEdge)
                        poolQ2(e).p(end + 1, 1) = pEdge;
                        poolQ2(e).map(end + 1, :) = [c g s t];
                    end
                catch
                end
            end

            % Store raw results; per-ERP BH-FDR post-pass fills q, sig, edges

            CompStats = storeCell(CompStats, erp, condName, tgtGrp, refGrp, IRR, CIlo, CIhi, pMat, pFDR_Mat, sig, emptyEdges(), NaN);
        end
    end
end

%% -------------------- Per-ERP BH-FDR post-pass -------------------- %%

for e = 1: numel(ERPs)
    erp = ERPs{e};

    pPool = poolQ2(e).p;
    mFDR = numel(pPool);

    if mFDR <= 0
        continue;
    end

    pFDR_Pool = bh_adjust(pPool);

    % Assign edge-level pFDR back into each panel matrix

    for ii = 1: mFDR
        c = poolQ2(e).map(ii, 1);
        g = poolQ2(e).map(ii, 2);
        s = poolQ2(e).map(ii, 3);
        t = poolQ2(e).map(ii, 4);

        condName = Conds{c};
        tgtGrp = char(TargetList{g});

        if isfield(CompStats, erp) && isfield(CompStats.(erp), condName) && isfield(CompStats.(erp).(condName), tgtGrp)
            CompStats.(erp).(condName).(tgtGrp).pFDR(s, t) = pFDR_Pool(ii);
            CompStats.(erp).(condName).(tgtGrp).mFDR = mFDR;
        end
    end

    % Build sig masks and edge lists per panel

    for c = 1: numel(Conds)
        condName = Conds{c};

        for g = 1: numel(TargetList)
            tgtGrp = char(TargetList{g});

            if ~isfield(CompStats, erp) || ~isfield(CompStats.(erp), condName) || ~isfield(CompStats.(erp).(condName), tgtGrp)
                continue;
            end

            IRR = CompStats.(erp).(condName).(tgtGrp).IRR;
            CIlo = CompStats.(erp).(condName).(tgtGrp).IRR_CIlo;
            CIhi = CompStats.(erp).(condName).(tgtGrp).IRR_CIhi;
            pMat = CompStats.(erp).(condName).(tgtGrp).p;
            pFDR_Mat = CompStats.(erp).(condName).(tgtGrp).pFDR;

            sig = isfinite(pFDR_Mat) & (pFDR_Mat <= 0.05);

            [sSig, tSig] = find(sig);

            clear source_dur target_dur

            for pop = 1: length(sSig)
                source_dur(pop) = nanmean(table2array(Table_target(Table_target.src == sSig(pop) & Table_target.tgt == tSig(pop), 8)));
                target_dur(pop) = nanmean(table2array(Table_target(Table_target.src == sSig(pop) & Table_target.tgt == tSig(pop), 9)));
            end

            if ~isempty(sSig)
                rr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                lo_vec = CIlo(sub2ind([7, 7], sSig, tSig));
                hi_vec = CIhi(sub2ind([7, 7], sSig, tSig));
                p_vec = pMat(sub2ind([7, 7], sSig, tSig));
                pFDR_vec = pFDR_Mat(sub2ind([7, 7], sSig, tSig));

                edges = table(sSig, tSig, rr_vec, lo_vec, hi_vec, p_vec, pFDR_vec, source_dur', target_dur', 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'p', 'pFDR', 'source_dur', 'target_dur'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);

                dlab = strings(numel(sSig), 1);
                dlab(rr_vec > 1) = "Higher vs " + string(refGrp);
                dlab(rr_vec < 1) = "Lower vs " + string(refGrp);
                edges.Direction = dlab;
            else
                edges = emptyEdges();
            end

            CompStats.(erp).(condName).(tgtGrp).sig = sig;
            CompStats.(erp).(condName).(tgtGrp).edges = edges;
            CompStats.(erp).(condName).(tgtGrp).mFDR = mFDR;
        end
    end
end

%% ------------------------------ Helpers ------------------------------ %%

    function tf = badSubject(S)
        tf = isempty(S) || ~isfield(S, 'trans_counts') || isempty(S.trans_counts) || S.M_segments <= 0 || S.N_transitions <= 0;
    end

    function [yvec, evec] = obs_vs_exp(S, srcIdx_local, tgtIdx_local)

        P = (S.C_segments(:) / S.M_segments);
        Ssq = sum(P .^ 2);
        den = 1 - Ssq;
        kK = numel(S.C_segments);

        if den <= eps
            yvec = NaN(size(srcIdx_local));
            evec = NaN(size(srcIdx_local));
            return;
        end

        gamma = 1 / den;

        E = gamma * S.N_transitions * (P * P.');
        E(1: kK + 1: kK * kK) = NaN;

        yvec = S.trans_counts(sub2ind([kK, kK], srcIdx_local, tgtIdx_local));
        evec = E(sub2ind([kK, kK], srcIdx_local, tgtIdx_local));
    end

    function S = getfield_safe(S, keys)

        for kk = 1: numel(keys)
            k = keys{kk};

            if ~isfield(S, k)
                S = [];
                return;
            end

            S = S.(k);

            if isempty(S)
                return;
            end
        end
    end

    function edges = emptyEdges()
        edges = table([], [], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'p', 'pFDR'});
        edges.srcLabel = strings(0, 1);
        edges.tgtLabel = strings(0, 1);
        edges.Direction = strings(0, 1);
    end

    function CompStats_local = storeCell(CompStats_local, erp_local, cond_local, tgtGrp_local, refGrp_local, IRR_local, CIlo_local, CIhi_local, pMat_local, pFDR_Mat_local, sig_local, edges_local, mFDR_local)

        if ~isfield(CompStats_local, erp_local)
            CompStats_local.(erp_local) = struct();
        end

        if ~isfield(CompStats_local.(erp_local), cond_local)
            CompStats_local.(erp_local).(cond_local) = struct();
        end

        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).IRR = IRR_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).IRR_CIlo = CIlo_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).IRR_CIhi = CIhi_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).p = pMat_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).pFDR = pFDR_Mat_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).sig = sig_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).edges = edges_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).RefGroup = refGrp_local;
        CompStats_local.(erp_local).(cond_local).(tgtGrp_local).mFDR = mFDR_local;
    end
end

function q = bh_adjust(pvals)

pvals = pvals(:);
[ps, ord] = sort(pvals, 'ascend');
m = numel(ps);
ranks = (1: m)';

qtmp = (m ./ ranks) .* ps;
qtmp = flipud(cummin(flipud(qtmp)));

q = NaN(size(ps));
q(ord) = min(1, qtmp);
end
