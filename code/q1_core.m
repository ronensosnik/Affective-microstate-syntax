function TSC = q1_core(Out, AgeMap, ERPs, Groups, Conds, ...
    FDR_Q1_SCOPE, GLM_OPTS, LOGE_CLAMP, MIN_OBS, MIN_EXP, MicrostateLabels, ...
    spec, alphaFDR, MedMap)
% spec: 'age_only' | 'classes' | 'poly'
[srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
TSC = struct();

for e = 1:numel(ERPs)
    erp = ERPs{e};
    for g = 1:numel(Groups)
        grp = Groups{g};

        % z-scored age within group
        Ages_raw = AgeMap.(grp)(:);
        muA = mean(Ages_raw, 'omitnan');
        sdA = std(Ages_raw, 'omitnan');
        if sdA <= eps, AgesZ = zeros(size(Ages_raw)); else, AgesZ = (Ages_raw - muA)./sdA; end

        % meds per subject (if applicable)
        if ~strcmp(spec,'age_only')
            switch spec
                case 'classes'
                    AD  = MedMap.(grp).AD(:); AP  = MedMap.(grp).AP(:);
                    MS  = MedMap.(grp).MS(:); ANX = MedMap.(grp).ANX(:);
                    OTH = MedMap.(grp).OTHER(:);
                case 'poly'
                    PC  = MedMap.(grp).poly_count(:);
            end
        end

        for c = 1:numel(Conds)
            condName  = Conds{c};
            condField = lower(condName);
            cellOut   = Out.(erp).(condField).(grp);

            IRR   = NaN(7);  beta0 = NaN(7);  se = NaN(7);  pval = NaN(7);
            ciLo  = NaN(7);  ciHi  = NaN(7);  nRows = zeros(7);  ObsOverExp = NaN(7);
            RowsAll = table(); % y, logE, covars..., src, tgt

            if ~isempty(cellOut)
                for sj = 1:numel(cellOut)
                    if sj > numel(AgesZ), continue; end
                    S = cellOut{sj};
                    if isempty(S) || ~isfield(S,'trans_counts') || isempty(S.trans_counts), continue; end
                    if S.M_segments <= 0 || S.N_transitions <= 0, continue; end

                    P = (S.C_segments(:) / S.M_segments);
                    E = S.N_transitions*(P * P.');
                    yvec = S.trans_counts(sub2ind([7,7], srcIdx, tgtIdx));
                    evec = E(sub2ind([7,7], srcIdx, tgtIdx));
                    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
                    if ~any(keep), continue; end

                    logE = log(evec(keep));
                    logE = max(min(logE, LOGE_CLAMP), -LOGE_CLAMP);

                    switch spec
                        case 'age_only'
                            R = table( yvec(keep), logE, repmat(AgesZ(sj), nnz(keep),1), ...
                                srcIdx(keep), tgtIdx(keep), ...
                                'VariableNames', {'y','logE','AgeZ','src','tgt'} );
                        case 'classes'
                            R = table( yvec(keep), logE, repmat(AgesZ(sj), nnz(keep),1), ...
                                repmat(AD(sj),  nnz(keep),1), repmat(AP(sj),  nnz(keep),1), ...
                                repmat(MS(sj),  nnz(keep),1), repmat(ANX(sj), nnz(keep),1), ...
                                repmat(OTH(sj), nnz(keep),1), ...
                                srcIdx(keep), tgtIdx(keep), ...
                                'VariableNames', {'y','logE','AgeZ','AD','AP','MS','ANX','OTHER','src','tgt'} );
                        case 'poly'
                            R = table( yvec(keep), logE, repmat(AgesZ(sj), nnz(keep),1), ...
                                repmat(PC(sj),  nnz(keep),1), ...
                                srcIdx(keep), tgtIdx(keep), ...
                                'VariableNames', {'y','logE','AgeZ','poly_count','src','tgt'} );
                    end
                    RowsAll = [RowsAll; R]; %#ok<AGROW>
                end
            end

            % Fit 42 edge GLMs
            for k = 1:numel(srcIdx)
                s = srcIdx(k); t = tgtIdx(k);
                R = RowsAll(RowsAll.src==s & RowsAll.tgt==t, :);
                if isempty(R), continue; end
                if sum(R.y) < MIN_OBS || sum(exp(R.logE)) < MIN_EXP, continue; end

                nRows(s,t) = height(R);
                ObsOverExp(s,t) = sum(R.y) / sum(exp(R.logE));

                try
                    switch spec
                        case 'age_only'
                            mdl = fitglm(R, 'y ~ 1 + AgeZ', ...
                                'Distribution','poisson','Link','log', ...
                                'Offset', R.logE, 'Options', GLM_OPTS);
                        case 'classes'
                            mdl = fitglm(R, 'y ~ 1 + AgeZ + AD + AP + MS + ANX + OTHER', ...
                                'Distribution','poisson','Link','log', ...
                                'Offset', R.logE, 'Options', GLM_OPTS);
                        case 'poly'
                            mdl = fitglm(R, 'y ~ 1 + AgeZ + poly_count', ...
                                'Distribution','poisson','Link','log', ...
                                'Offset', R.logE, 'Options', GLM_OPTS);
                    end
                    b   = mdl.Coefficients.Estimate;
                    seB = mdl.Coefficients.SE;
                    pB  = mdl.Coefficients.pValue;
                    CI  = coefCI(mdl, 0.05);

                    beta0(s,t) = b(1);
                    se(s,t)    = seB(1);
                    pval(s,t)  = pB(1);
                    IRR(s,t)   = exp(b(1));
                    ciLo(s,t)  = exp(CI(1,1));
                    ciHi(s,t)  = exp(CI(1,2));
                catch
                    % leave NaNs on failure
                end
            end

            % FDR within this cell (42 off-diagonal)
            if strcmp(FDR_Q1_SCOPE,'per_cell_42')
                [sigMask, qMat] = bh_fdr_mask(pval, alphaFDR);
            else
                sigMask = false(7); qMat = NaN(7);
            end

            dirMatrix = zeros(7);
            dirMatrix(sigMask & IRR > 1) =  1;
            dirMatrix(sigMask & IRR < 1) = -1;

            % Edge list (for convenience)
            [sSig, tSig] = find(sigMask);
            if ~isempty(sSig)
                irr_vec = IRR(sub2ind([7,7], sSig, tSig));
                lo_vec  = ciLo(sub2ind([7,7], sSig, tSig));
                hi_vec  = ciHi(sub2ind([7,7], sSig, tSig));
                q_vec   = qMat(sub2ind([7,7], sSig, tSig));
                o2e_vec = ObsOverExp(sub2ind([7,7], sSig, tSig));
                edges = table(sSig, tSig, irr_vec, lo_vec, hi_vec, o2e_vec, q_vec, ...
                    'VariableNames',{'src','tgt','IRR','IRR_CIlo','IRR_CIhi','ObsOverExp','qFDR'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [],1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [],1);
                lab = strings(numel(sSig),1);
                lab(irr_vec>1)="Above expected"; lab(irr_vec<1)="Below expected";
                edges.Direction = lab;
            else
                edges = table([],[],[],[],[],[],[], ...
                    'VariableNames',{'src','tgt','IRR','IRR_CIlo','IRR_CIhi','ObsOverExp','qFDR'});
                edges.srcLabel = strings(0,1);
                edges.tgtLabel = strings(0,1);
                edges.Direction = strings(0,1);
            end

            TSC.(erp).(grp).(condName).IRR = IRR;
            TSC.(erp).(grp).(condName).beta0 = beta0;
            TSC.(erp).(grp).(condName).se = se;
            TSC.(erp).(grp).(condName).p = pval;
            TSC.(erp).(grp).(condName).q = qMat;
            TSC.(erp).(grp).(condName).sig = sigMask;
            TSC.(erp).(grp).(condName).dir = dirMatrix;
            TSC.(erp).(grp).(condName).nRows = nRows;
            TSC.(erp).(grp).(condName).ObsOverExp = ObsOverExp;
            TSC.(erp).(grp).(condName).FDR_SCOPE = FDR_Q1_SCOPE;
            TSC.(erp).(grp).(condName).edges = edges;
        end
    end
end
end