function [rt_results, rt_simple_slopes] = fitConfirmatoryRT(Tall, mode_str, with_interaction)
% fitConfirmatoryRT:  Confirmatory models for RATE or DELTA predictors.
%
% mode_str:  'rate' or 'delta'
% with_interaction:  true to include Predictor×Condition and compute simple slopes
%

rows = {};
rowsSS = {};

erp_list = {'N200', 'P300', 'LPP'};

for e = 1: numel(erp_list)
    erp = erp_list{e};
    T = Tall(Tall.ERP == string(erp), :);

    if strcmp(erp, 'LPP')
        base = 'backbone_';
        preds = {[base mode_str]};
        pnames = {'LPP_backbone'};
    else
        base1 = 'anchor1_';
        base2 = 'anchor2_';
        preds = {[base1 mode_str], [base2 mode_str]};
        if strcmp(erp, 'N200')
            pnames = {'N200_anchor1', 'N200_anchor2'};
        else
            pnames = {'P300_anchor1', 'P300_anchor2'};
        end
    end

    for p = 1: numel(preds)
        pr = preds{p};
        pname = pnames{p};

        zcol = ['z_' pr];

        mu = mean(T.(pr), 'omitnan');
        sd = std(T.(pr), 'omitnan');

        if sd == 0 || isnan(sd)
            T.(zcol) = zeros(height(T), 1);
        else
            T.(zcol) = (T.(pr) - mu) ./ sd;
        end

        if with_interaction
            formula = ['log_median_rt ~ ' zcol ' + condition + ' zcol ':condition + age + (1|subject)'];
        else
            formula = ['log_median_rt ~ ' zcol ' + condition + age + (1|subject)'];
        end

        use_rows = isfinite(T.(zcol)) & isfinite(T.log_median_rt) & isfinite(double(T.age));

        TM = T(use_rows, :);

        lme = fitlme(TM, formula, 'FitMethod', 'REML');

        coef = lme.Coefficients;
        row_idx = strcmp(string(coef.Name), string(zcol));

        beta = coef.Estimate(row_idx);
        ci_low = coef.Lower(row_idx);
        ci_high = coef.Upper(row_idx);
        pval = coef.pValue(row_idx);

        r = table;
        r.ERP = string(erp);
        r.Predictor = string([pname '_' mode_str]);
        r.BetaStd = beta;
        r.CI_low = ci_low;
        r.CI_high = ci_high;
        r.pValue = pval;

        if with_interaction
            % Overall interaction p-value (all interaction terms jointly)
            term_names = string(coef.Name);
            mask_int = startsWith(term_names, string([zcol ':condition_']));
            if any(mask_int)
                H = zeros(sum(mask_int), height(coef));
                cols = find(mask_int);
                for i = 1: numel(cols)
                    H(i, cols(i)) = 1;
                end
                p_inter = coefTest(lme, H, zeros(sum(mask_int), 1));
            else
                p_inter = NaN;
            end
            r.pValue_interaction = p_inter;

            % Simple slopes per condition using refit with that condition as reference
            cond_levels = categories(TM.condition);
            for i = 1: numel(cond_levels)
                cond = string(cond_levels{i});
                [b, ci, p] = simpleSlopeRefit(TM, formula, zcol, cond);
                rs = table;
                rs.ERP = string(erp);
                rs.Predictor = string([pname '_' mode_str]);
                rs.Condition = cond;
                rs.BetaStd = b;
                rs.CI_low = ci(1);
                rs.CI_high = ci(2);
                rs.pValue = p;
                rowsSS{end + 1, 1} = rs;
            end
        else
            r.pValue_interaction = NaN;
        end

        rows{end + 1, 1} = r;
    end
end

rt_results = vertcat(rows{ : });

rt_results.qValue = fdr_bh_local(rt_results.pValue);

if with_interaction
    rt_simple_slopes = vertcat(rowsSS{ : });
    rt_simple_slopes.qValue = fdr_bh_local(rt_simple_slopes.pValue);
else
    rt_simple_slopes = table;
end

end
