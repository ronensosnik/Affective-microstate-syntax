function [out_main, out_slopes] = fitTrialRT_Interactions(T)

% fitTrialRT_Interactions: Trial-level LME with Predictor×Condition and random slopes by subject.
%
% Inputs:
%   T: table with columns
%      subject, condition (Neutral/Negative/Positive, categorical), age, trial_z,
%      z_N200_a1_rate ... z_LPP_backbone_rate, log_rt
%
% Outputs:
%   out_main   : one row per predictor (slope @ Neutral + joint interaction p; BH across predictors)
%   out_slopes : simple slopes per condition (BH across all simple slopes)

% Make sure condition order is Neutral, Negative, Positive

T.condition = categorical(T.condition, {'Neutral', 'Negative', 'Positive'});

% Keep complete cases for all variables we will touch

preds = {'z_N200_a1_rate', 'z_N200_a2_rate', 'z_P300_a1_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
use = isfinite(T.log_rt) & isfinite(double(T.age)) & isfinite(T.trial_z);

for k = 1: numel(preds)
    use = use & isfinite(T.(preds{k}));
end

T = T(use, :);

% Fixed-effect formula (main effects + Predictor×Condition for all five)

fix = ['log_rt ~ z_N200_a1_rate*condition + z_N200_a2_rate*condition + ' 'z_P300_a1_rate*condition + z_P300_a2_rate*condition + z_LPP_backbone_rate*condition + ' 'age + trial_z' ];

% Random-effects (try rich model first; then fall back if needed)

rand_full = '(1 + z_N200_a1_rate + z_N200_a2_rate + z_P300_a1_rate + z_P300_a2_rate + z_LPP_backbone_rate | subject)';
rand_mid = '(1 + z_N200_a2_rate + z_P300_a2_rate | subject)';
rand_min = '(1 | subject)';

formula_try = [fix ' + ' rand_full];

try
    lme = fitlme(T, formula_try, 'FitMethod', 'REML', 'CovariancePattern', 'Diagonal');
catch
    try
        formula_try = [fix ' + ' rand_mid];
        lme = fitlme(T, formula_try, 'FitMethod', 'REML', 'CovariancePattern', 'Diagonal');
    catch
        formula_try = [fix ' + ' rand_min];
        lme = fitlme(T, formula_try, 'FitMethod', 'REML');
    end
end

C = lme.Coefficients;
names = string(C.Name);

% Slope @ Neutral (reference) for each predictor

want = {'z_N200_a1_rate', 'z_N200_a2_rate', 'z_P300_a1_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
pretty = {'N200 anchor1', 'N200 anchor2', 'P300 anchor1', 'P300 anchor2', 'LPP backbone'};

rows = {};
p_int = nan(numel(want), 1);

for i = 1: numel(want)
    zcol = string(want{i});
    r = table;

    r.Predictor = string(pretty{i});

    idx = strcmp(names, zcol);
    r.BetaStd_Neutral = C.Estimate(idx);
    r.CI_low = C.Lower(idx);
    r.CI_high = C.Upper(idx);
    r.pValue = C.pValue(idx);

    % Joint interaction test for this predictor (Negative and Positive vs Neutral)

    termN = string(zcol + ':condition_Negative');
    termP = string(zcol + ':condition_Positive');

    mask = ismember(names, [termN termP]);

    if any(mask)
        H = zeros(sum(mask), height(C));
        cols = find(mask);

        for j = 1: numel(cols)
            H(j, cols(j)) = 1;
        end
        p_inter = coefTest(lme, H, zeros(sum(mask), 1));
    else
        p_inter = NaN;
    end

    r.p_interaction = p_inter;
    p_int(i) = p_inter;

    rows{end + 1, 1} = r;
end

out_main = vertcat(rows{ : });

% BH–FDR across the five Neutral slopes

out_main.qValue = fdr_bh_local(out_main.pValue);

% -------- Simple slopes by condition via refit with that condition as reference

rowsS = {};
conds = {'Neutral','Negative','Positive'};

for i = 1: numel(want)
    zcol = string(want{i});
    lab  = string(pretty{i});

    for c = 1: numel(conds)
        cond = string(conds{c});
        [b, ci, p] = simpleSlopeRefit(T, [fix ' + ' rand_min], zcol, cond); % refit with random intercept for stability
        rs = table;
        rs.Predictor = lab;
        rs.Condition = cond;
        rs.BetaStd = b;
        rs.CI_low = ci(1);
        rs.CI_high = ci(2);
        rs.pValue = p;
        rowsS{end + 1, 1} = rs;
    end
end

out_slopes = vertcat(rowsS{ : });

% BH–FDR across all simple slopes (5 predictors × 3 conditions = 15 tests)

out_slopes.qValue = fdr_bh_local(out_slopes.pValue);

end