function CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, varargin)

% buildAnchoredContrasts (per-target 42-FDR version)
% Stage 1 (Ref-only, per pair): y ~ Poisson, log(mu) = log(E) + α0 + α_age*AgeS_ref
% Stage 2 (Ref+Target, per pair): y ~ Poisson, log(mu) = [log(E)+α0+α_age*AgeS_ref(row)] + δ*G (no intercept)
% FDR: BH within each ERP×Condition×TargetGroup across the 42 pairs.

p = inputParser;
addParameter(p, 'RefGroup', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'TargetGroups', {}, @(x)iscell(x)||isstring(x));
addParameter(p, 'Alpha', 0.05);
addParameter(p, 'MaxIter', 2000);
addParameter(p, 'MinObs', 3);
addParameter(p, 'MinExp', 3);
addParameter(p, 'LogEClamp', 12);
parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.RefGroup), 'Please provide ''RefGroup'', e.g., ''HC''.');
assert(any(strcmpi(Groups, opt.RefGroup)), 'RefGroup "%s" not found in Groups.', char(opt.RefGroup));

% Target list

if isempty(opt.TargetGroups)
    TargetList = Groups(~strcmpi(Groups, opt.RefGroup)); % all except reference
else
    tg = cellstr(opt.TargetGroups);
    TargetList = tg(~strcmpi(tg, opt.RefGroup) & ismember(tg, Groups));
end

assert(~isempty(TargetList), 'No valid TargetGroups to compare against RefGroup.');

[srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
GLM_OPTS = statset('MaxIter', opt.MaxIter, 'TolX', 1e-10, 'TolFun', 1e-10);
CompStats = struct();
MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

for e = 1: numel(ERPs)
    erp = ERPs{e};

    for c = 1: numel(Conds)
        condName = Conds{c};
        condKey = lower(condName);

        RefCell = getfield_safe(Out, {erp, condKey, char(opt.RefGroup)});

        % Prebuild ref rows ONCE (shared across targets)

        Tref = table();
        if ~isempty(RefCell)
            agesRef = AgeMap.(char(opt.RefGroup));

            for sj = 1: numel(RefCell)
                S = RefCell{sj};

                if badSubject(S)
                    continue;
                end

                [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                if ~any(keep)
                    continue;
                end

                logE = log(evec(keep));
                logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                R = table(yvec(keep), logE, repmat(agesRef(min(sj, end)), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeRaw', 'src', 'tgt'});
                Tref = [Tref; R];
            end
        end

        for g = 1: numel(TargetList)
            tgtGrp = TargetList{g};
            Ttgt = table();
            TgtCell = getfield_safe(Out, {erp, condKey, tgtGrp});

            if ~isempty(TgtCell)
                agesTgt = AgeMap.(tgtGrp);

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

                    logE = log(evec(keep));
                    logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                    R = table(yvec(keep), logE, repmat(agesTgt(min(sj, end)), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeRaw', 'src', 'tgt'});
                    Ttgt = [Ttgt; R];
                end
            end

            % Init per-target containers

            IRR = NaN(7);
            CIlo = NaN(7);
            CIhi = NaN(7);
            pMat = NaN(7);
            qMat = NaN(7);
            sig = false(7);

            if isempty(Tref) || isempty(Ttgt)
                store(erp, condName, tgtGrp, opt.RefGroup, IRR, CIlo, CIhi, pMat, qMat, sig, emptyEdges());
                continue;
            end

            % --- Fit per pair

            for k = 1: numel(srcIdx)
                s = srcIdx(k);
                t = tgtIdx(k);

                Rref = Tref(Tref.src == s & Tref.tgt == t, :);
                Rtgt = Ttgt(Ttgt.src == s & Ttgt.tgt == t, :);

                if isempty(Rref) || isempty(Rtgt)
                    continue;
                end

                % Stage 1: reference-only fit (age z-scored within ref)

                muR = mean(Rref.AgeRaw, 'omitnan');
                sdR = std(Rref.AgeRaw, 'omitnan');
                Rref_fit = Rref;

                if sdR <= eps
                    Rref_fit.AgeS = zeros(height(Rref_fit), 1);
                else
                    Rref_fit.AgeS = (Rref_fit.AgeRaw - muR)./sdR;
                end

                % screens for ref fit

                if sum(Rref_fit.y) < opt.MinObs || sum(exp(Rref_fit.logE)) < opt.MinExp
                    continue;
                end

                try
                    mdlRef = fitglm(Rref_fit, 'y ~ 1 + AgeS', 'Distribution', 'poisson', 'Link', 'log', 'Offset', Rref_fit.logE, 'Options', GLM_OPTS);

                    alpha0 = mdlRef.Coefficients.Estimate(1);
                    alphaAge = 0;

                    if any(strcmp(mdlRef.CoefficientNames, 'AgeS'))
                        alphaAge = mdlRef.Coefficients.Estimate(strcmp(mdlRef.CoefficientNames, 'AgeS'));
                    end
                catch
                    continue;
                end

                % Stage 2: combine ref+target & anchor offset at each row's age

                Rref2 = dropIfExists(Rref, {'AgeS', 'AgeS_ref', 'offset_ref', 'Gnum'});
                Rtgt2 = dropIfExists(Rtgt, {'AgeS', 'AgeS_ref', 'offset_ref', 'Gnum'});
                R2 = [Rref2; Rtgt2];

                if sdR <= eps
                    R2.AgeS_ref = zeros(height(R2), 1);
                else
                    R2.AgeS_ref = (R2.AgeRaw - muR)./sdR;
                end

                R2.offset_ref = R2.logE + alpha0 + alphaAge * R2.AgeS_ref;
                R2.Gnum = [zeros(height(Rref2), 1); ones(height(Rtgt2), 1)];

                % screens for combined fit

                if sum(R2.y) < opt.MinObs || sum(exp(R2.offset_ref)) < opt.MinExp
                    continue;
                end

                % quasi-separation guard

                y_t = sum(R2.y(R2.Gnum == 1));
                y_r = sum(R2.y(R2.Gnum == 0));

                if (y_t == 0 || y_r == 0)
                    continue;
                end

                try
                    mdl2 = fitglm(R2, 'y ~ -1 + Gnum', 'Distribution', 'poisson', 'Link', 'log', 'Offset', R2.offset_ref, 'Options', GLM_OPTS);

                    b = mdl2.Coefficients.Estimate(1);
                    p = mdl2.Coefficients.pValue(1);
                    CI = coefCI(mdl2, 0.05);

                    IRR(s, t)  = exp(b);
                    CIlo(s, t) = exp(CI(1, 1));
                    CIhi(s, t) = exp(CI(1, 2));
                    pMat(s, t) = p;

                catch
                    % leave NaNs
                end
            end

            % --- BH–FDR within this Target × ERP × Condition across 42 pairs

            offdiag = ~eye(7);
            p_lin = pMat(offdiag);
            keep = ~isnan(p_lin);

            if any(keep)
                q_all = bh_adjust(p_lin(keep));
                q_only = NaN(size(p_lin));
                q_only(keep) = q_all;

                qMat(offdiag) = q_only;
            end

            sig = qMat < opt.Alpha & ~isnan(qMat);

            % --- Edge list

            [sSig, tSig] = find(sig);

            if ~isempty(sSig)
                rr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                lo_vec = CIlo(sub2ind([7, 7], sSig, tSig));
                hi_vec = CIhi(sub2ind([7, 7], sSig, tSig));
                q_vec  = qMat(sub2ind([7, 7], sSig, tSig));

                edges = table(sSig, tSig, rr_vec, lo_vec, hi_vec, q_vec, 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);
                dlab = strings(numel(sSig), 1);
                dlab(rr_vec > 1) = "Higher vs " + string(opt.RefGroup);
                dlab(rr_vec < 1) = "Lower vs " + string(opt.RefGroup);
                edges.Direction = dlab;
            else
                edges = emptyEdges();
            end

            % --- Store

            store(erp, condName, tgtGrp, opt.RefGroup, IRR, CIlo, CIhi, pMat, qMat, sig, edges);
        end

        CompStats.(erp).(condName).FDR_SCOPE = 'per_group_condition_42pairs';
    end
end

    function tf = badSubject(S)
        tf = isempty(S) || ~isfield(S, 'trans_counts') || isempty(S.trans_counts) || S.M_segments <= 0 || S.N_transitions <= 0;
    end

    function [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx)
        P = (S.C_segments(:) / S.M_segments); % 7×1
        E = S.N_transitions * (P * P.'); % 7×7 under independence
        yvec = S.trans_counts(sub2ind([7, 7], srcIdx, tgtIdx));
        evec = E(sub2ind([7, 7], srcIdx, tgtIdx));
    end

    function T = dropIfExists(T, names)
        keep = ismember(names, T.Properties.VariableNames);

        if any(keep)
            T = removevars(T, names(keep));
        end
    end

    function edges = emptyEdges()
        edges = table([], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});
        edges.srcLabel = strings(0, 1);
        edges.tgtLabel = strings(0, 1);
        edges.Direction = strings(0, 1);
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

    function store(erp, condName, tgtGrp, refGrp, IRR, CIlo, CIhi, pMat, qMat, sig, edges)

        persistent CompStats_local

        if isempty(CompStats_local)
            CompStats_local = struct();
        end

        CompStats_local.(erp).(condName).(tgtGrp).IRR = IRR;
        CompStats_local.(erp).(condName).(tgtGrp).IRR_CIlo = CIlo;
        CompStats_local.(erp).(condName).(tgtGrp).IRR_CIhi = CIhi;
        CompStats_local.(erp).(condName).(tgtGrp).p = pMat;
        CompStats_local.(erp).(condName).(tgtGrp).q = qMat;
        CompStats_local.(erp).(condName).(tgtGrp).sig = sig;
        CompStats_local.(erp).(condName).(tgtGrp).edges = edges;
        CompStats_local.(erp).(condName).(tgtGrp).RefGroup = char(refGrp);
        assignin('caller', 'CompStats', CompStats_local);
    end

    function q = bh_adjust(p)
        p = p(:);
        [ps, ord] = sort(p, 'ascend');
        m = numel(ps);
        ranks = (1: m)';
        qtmp = (m ./ ranks) .* ps;
        qtmp = flipud(cummin(flipud(qtmp)));
        q = NaN(size(ps));
        q(ord) = min(1, qtmp);
    end
end
