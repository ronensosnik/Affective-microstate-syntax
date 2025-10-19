function edge_results = fitEdgeRT(Tall)
% fitEdgeRT:  Mixed-effects per edge with BH FDR within ERP.

erp_list = {'N200', 'P300', 'LPP'};

rows = {};

for e = 1: numel(erp_list)
    erp = erp_list{e};
    TE = Tall(Tall.ERP == string(erp), :);

    % z-standardize per ERP
    mu = mean(TE.edge_rate, 'omitnan');
    sd = std(TE.edge_rate, 'omitnan');
    if sd == 0 || isnan(sd)
        TE.z_edge = zeros(height(TE), 1);
    else
        TE.z_edge = (TE.edge_rate - mu) ./ sd;
    end

    edges = unique(TE.Edge);

    pvals = zeros(numel(edges), 1);

    betas = zeros(numel(edges), 1);

    ciL = zeros(numel(edges), 1);

    ciH = zeros(numel(edges), 1);

    for k = 1: numel(edges)
        ed = edges(k);
        TK = TE(TE.Edge == ed, :);
        use_rows = isfinite(TK.z_edge) & isfinite(TK.log_median_rt) & isfinite(double(TK.age));
        TK = TK(use_rows, :);

        if height(TK) < 20
            betas(k) = NaN;
            ciL(k) = NaN;
            ciH(k) = NaN;
            pvals(k) = NaN;
            continue
        end

        formula = 'log_median_rt ~ z_edge + condition + age + (1|subject)';

        lme = fitlme(TK, formula, 'FitMethod', 'REML');

        coef = lme.Coefficients;
        row_idx = strcmp(string(coef.Name), 'z_edge');
        betas(k) = coef.Estimate(row_idx);
        ciL(k) = coef.Lower(row_idx);
        ciH(k) = coef.Upper(row_idx);
        pvals(k) = coef.pValue(row_idx);
    end

    q = fdr_bh_local(pvals);

    Tout = table;
    Tout.ERP = repmat(string(erp), numel(edges), 1);
    Tout.Edge = edges;
    Tout.BetaStd = betas;
    Tout.CI_low = ciL;
    Tout.CI_high = ciH;
    Tout.pValue = pvals;
    Tout.qValue = q;

    rows{end + 1, 1} = Tout;
end

edge_results = vertcat(rows{ : });

end
