function CompStats = buildAnchoredContrasts_meds(Out, AgeMap, MedMap, ERPs, Conds, Groups, varargin)
% Anchored contrasts with medication covariates in Stage 2.
% Stage 1 (Reference-only):   y ~ Poisson, log(mu) = log(E) + α0 + α_age*AgeS_ref
% Stage 2 (Ref+Target):       y ~ Poisson, log(mu) = [log(E)+α0+α_age*AgeS_ref(row)] + δ*G + meds
% FDR: BH within each ERP×Condition×Target across the 42 edges.

p = inputParser;
addParameter(p, 'RefGroup', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'TargetGroups', {}, @(x)iscell(x)||isstring(x));
addParameter(p, 'Spec', 'classes', @(s) any(strcmp(s,{'classes','poly'})));
addParameter(p, 'Alpha', 0.05);
addParameter(p, 'MaxIter', 2000);
addParameter(p, 'MinObs', 3);
addParameter(p, 'MinExp', 3);
addParameter(p, 'LogEClamp', 12);
parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.RefGroup), 'Please provide ''RefGroup'', e.g., ''HC''.');
assert(any(strcmpi(Groups, opt.RefGroup)), 'RefGroup "%s" not found in Groups.', char(opt.RefGroup));

if isempty(opt.TargetGroups)
    TargetList = Groups(~strcmpi(Groups, opt.RefGroup));
else
    tg = cellstr(opt.TargetGroups);
    TargetList = tg(~strcmpi(tg, opt.RefGroup) & ismember(tg, Groups));
end
assert(~isempty(TargetList), 'No valid TargetGroups to compare against RefGroup.');

[srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
GLM_OPTS = statset('MaxIter', opt.MaxIter, 'TolX', 1e-10, 'TolFun', 1e-10);
CompStats = struct();
MicrostateLabels = {'A','B','C','D','E','F','G'};



for e = 1:numel(ERPs)
    erp = ERPs{e};
    for c = 1:numel(Conds)
        condName = Conds{c};
        condKey  = lower(condName);

        RefCell = nested_get(Out, {erp, condKey, char(opt.RefGroup)});
        Tref = table();

        if ~isempty(RefCell)
            agesRef = AgeMap.(char(opt.RefGroup));
            medsRef = MedMap.(char(opt.RefGroup));
            for sj = 1:numel(RefCell)
                S = RefCell{sj};
                if badSubject(S), continue; end
                [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                keep = isfinite(evec) & evec>0 & ~isnan(yvec);
                if ~any(keep), continue; end

                logE = log(evec(keep));
                logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                switch opt.Spec
                    case 'classes'
                        R = table(yvec(keep), logE, repmat(agesRef(min(sj,end)), nnz(keep),1), ...
                            repmat(medsRef.AD(min(sj,end)),    nnz(keep),1), ...
                            repmat(medsRef.AP(min(sj,end)),    nnz(keep),1), ...
                            repmat(medsRef.MS(min(sj,end)),    nnz(keep),1), ...
                            repmat(medsRef.ANX(min(sj,end)),   nnz(keep),1), ...
                            repmat(medsRef.OTHER(min(sj,end)), nnz(keep),1), ...
                            srcIdx(keep), tgtIdx(keep), ...
                            'VariableNames', {'y','logE','AgeRaw','AD','AP','MS','ANX','OTHER','src','tgt'});
                    case 'poly'
                        R = table(yvec(keep), logE, repmat(agesRef(min(sj,end)), nnz(keep),1), ...
                            repmat(medsRef.poly_count(min(sj,end)), nnz(keep),1), ...
                            srcIdx(keep), tgtIdx(keep), ...
                            'VariableNames', {'y','logE','AgeRaw','poly_count','src','tgt'});
                end
                Tref = [Tref; R]; %#ok<AGROW>
            end
        end

        for g = 1:numel(TargetList)
            tgtGrp = TargetList{g};
            Ttgt = table();
            TgtCell = nested_get(Out, {erp, condKey, tgtGrp});
            if ~isempty(TgtCell)
                agesT = AgeMap.(tgtGrp);
                medsT = MedMap.(tgtGrp);
                for sj = 1:numel(TgtCell)
                    S = TgtCell{sj};
                    if badSubject(S), continue; end
                    [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                    keep = isfinite(evec) & evec>0 & ~isnan(yvec);
                    if ~any(keep), continue; end

                    logE = log(evec(keep));
                    logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                    switch opt.Spec
                        case 'classes'
                            R = table(yvec(keep), logE, repmat(agesT(min(sj,end)), nnz(keep),1), ...
                                repmat(medsT.AD(min(sj,end)),    nnz(keep),1), ...
                                repmat(medsT.AP(min(sj,end)),    nnz(keep),1), ...
                                repmat(medsT.MS(min(sj,end)),    nnz(keep),1), ...
                                repmat(medsT.ANX(min(sj,end)),   nnz(keep),1), ...
                                repmat(medsT.OTHER(min(sj,end)), nnz(keep),1), ...
                                srcIdx(keep), tgtIdx(keep), ...
                                'VariableNames', {'y','logE','AgeRaw','AD','AP','MS','ANX','OTHER','src','tgt'});
                        case 'poly'
                            R = table(yvec(keep), logE, repmat(agesT(min(sj,end)), nnz(keep),1), ...
                                repmat(medsT.poly_count(min(sj,end)), nnz(keep),1), ...
                                srcIdx(keep), tgtIdx(keep), ...
                                'VariableNames', {'y','logE','AgeRaw','poly_count','src','tgt'});
                    end
                    Ttgt = [Ttgt; R]; %#ok<AGROW>
                end
            end

            IRR = NaN(7); CIlo = NaN(7); CIhi = NaN(7); pMat = NaN(7); qMat = NaN(7); sig = false(7);

            if isempty(Tref) || isempty(Ttgt)
                store_node(); continue;
            end

            % Stage 1: fit reference-only (age standardized within ref)
            muR = mean(Tref.AgeRaw, 'omitnan'); sdR = std(Tref.AgeRaw, 'omitnan');
            Rref = Tref;
            if sdR <= eps, Rref.AgeS = zeros(height(Rref),1);
            else,          Rref.AgeS = (Rref.AgeRaw - muR)./sdR;
            end

            for k = 1:numel(srcIdx)
                s = srcIdx(k); t = tgtIdx(k);

                Rref_k = Rref(Rref.src==s & Rref.tgt==t, :);
                if isempty(Rref_k), continue; end
                if sum(Rref_k.y) < opt.MinObs || sum(exp(Rref_k.logE)) < opt.MinExp, continue; end

                try
                    mdlRef = fitglm(Rref_k, 'y ~ 1 + AgeS', ...
                        'Distribution','poisson','Link','log', ...
                        'Offset', Rref_k.logE, 'Options', GLM_OPTS);
                    a0 = mdlRef.Coefficients.Estimate(1);
                    aA = 0; if any(strcmp(mdlRef.CoefficientNames,'AgeS'))
                        aA = mdlRef.Coefficients.Estimate(strcmp(mdlRef.CoefficientNames,'AgeS'));
                    end
                catch
                    continue;
                end

                % Stage 2: ref+target; anchor with ref offset at each row's age
                R2_ref = Tref(Tref.src==s & Tref.tgt==t, :);
                R2_tgt = Ttgt(Ttgt.src==s & Ttgt.tgt==t, :);
                if isempty(R2_ref) || isempty(R2_tgt), continue; end

                R2 = [R2_ref; R2_tgt];
                if sdR <= eps, R2.AgeS_ref = zeros(height(R2),1);
                else,          R2.AgeS_ref = (R2.AgeRaw - muR)./sdR;
                end
                R2.offset_ref = R2.logE + a0 + aA * R2.AgeS_ref;
                R2.Gnum       = [zeros(height(R2_ref),1); ones(height(R2_tgt),1)];

                % quasi-separation guard & screens
                if sum(R2.y) < opt.MinObs || sum(exp(R2.offset_ref)) < opt.MinExp, continue; end
                if sum(R2.y(R2.Gnum==0))==0 || sum(R2.y(R2.Gnum==1))==0, continue; end

                try
                    switch opt.Spec
                        case 'classes'
                            mdl2 = fitglm(R2, 'y ~ -1 + Gnum + AD + AP + MS + ANX + OTHER', ...
                                'Distribution','poisson','Link','log', ...
                                'Offset', R2.offset_ref, 'Options', GLM_OPTS);
                        case 'poly'
                            mdl2 = fitglm(R2, 'y ~ -1 + Gnum + poly_count', ...
                                'Distribution','poisson','Link','log', ...
                                'Offset', R2.offset_ref, 'Options', GLM_OPTS);
                    end
                    b  = mdl2.Coefficients.Estimate(strcmp(mdl2.CoefficientNames,'Gnum'));
                    p  = mdl2.Coefficients.pValue(strcmp(mdl2.CoefficientNames,'Gnum'));
                    CI = coefCI(mdl2, 0.05);

                    % When there are multiple coefficients, ensure we picked Gnum row:
                    if numel(b)~=1
                        % safest: locate exact name
                        gi = find(strcmp(mdl2.CoefficientNames,'Gnum'),1);
                        b  = mdl2.Coefficients.Estimate(gi);
                        p  = mdl2.Coefficients.pValue(gi);
                        CI = coefCI(mdl2,0.05);
                        CI = CI(gi,:);
                    end

                    IRR(s,t)  = exp(b);
                    CIlo(s,t) = exp(CI(1));
                    CIhi(s,t) = exp(CI(2));
                    pMat(s,t) = p;

                catch
                    % leave NaNs
                end
            end

            % BH–FDR within Target×ERP×Condition (42 pairs)
            offdiag = ~eye(7);
            p_lin = pMat(offdiag);
            keep = ~isnan(p_lin);
            if any(keep)
                q_all = bh_adjust(p_lin(keep));
                q_only = NaN(size(p_lin)); q_only(keep) = q_all;
                qMat(offdiag) = q_only;
            end
            sig = qMat < opt.Alpha & ~isnan(qMat);

            % Edge list (for convenience)
            [sSig, tSig] = find(sig);
            if ~isempty(sSig)
                rr = IRR(sub2ind([7,7], sSig, tSig));
                lo = CIlo(sub2ind([7,7], sSig, tSig));
                hi = CIhi(sub2ind([7,7], sSig, tSig));
                qv = qMat(sub2ind([7,7], sSig, tSig));
                edges = table(sSig, tSig, rr, lo, hi, qv, ...
                    'VariableNames',{'src','tgt','IRR','IRR_CIlo','IRR_CIhi','qFDR'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)),[],1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)),[],1);
                dlab = strings(numel(sSig),1);
                dlab(rr>1) = "Higher vs " + string(opt.RefGroup);
                dlab(rr<1) = "Lower vs "  + string(opt.RefGroup);
                edges.Direction = dlab;
            else
                edges = emptyEdges();
            end

            CompStats.(erp).(condName).(tgtGrp).IRR      = IRR;
            CompStats.(erp).(condName).(tgtGrp).IRR_CIlo  = CIlo;
            CompStats.(erp).(condName).(tgtGrp).IRR_CIhi  = CIhi;
            CompStats.(erp).(condName).(tgtGrp).p         = pMat;
            CompStats.(erp).(condName).(tgtGrp).q         = qMat;
            CompStats.(erp).(condName).(tgtGrp).sig       = sig;
            CompStats.(erp).(condName).(tgtGrp).edges     = edges;
            CompStats.(erp).(condName).(tgtGrp).RefGroup  = char(opt.RefGroup);
            CompStats.(erp).(condName).FDR_SCOPE = 'per_group_condition_42pairs';

        end
    end
end

end