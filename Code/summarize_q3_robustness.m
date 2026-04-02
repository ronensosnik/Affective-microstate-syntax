function OUT = summarize_q3_robustness(base_csv_or_table, sens_csv_or_table, varargin)

% SUMMARIZE_Q3_ROBUSTNESS
% Compare Q3 trial-level RT~syntax results between baseline and a sensitivity pipeline.
%
% Inputs:
%   base_csv_or_table: path to RT_trial_confirmatory_base.csv OR a table
%   sens_csv_or_table: path to RT_trial_confirmatory_*.csv OR a table
%
% Name-value options:
%   'BaseLabel'        (default: "baseline")
%   'SensLabel'        (default: "sensitivity")
%   'Alpha'            (default: 0.05)
%   'WriteCSV'         (default: true)  : write *_detail.csv and *_overall.csv
%   'WriteMAT'         (default: true)  : write *.mat containing OUT
%   'OutDir'           (default: inferred from sens CSV folder, else pwd)
%   'OutStem'          (default: inferred as "Q3_robustness_<Base>_vs_<Sens>")
%   'OutMAT'           (default: fullfile(OutDir, OutStem + ".mat"))
%   'OutDetailCSV'     (default: fullfile(OutDir, OutStem + "_detail.csv"))
%   'OutOverallCSV'    (default: fullfile(OutDir, OutStem + "_overall.csv"))
%
% Output struct OUT:
%   OUT.detail   : table (N200 Anchor-2, P300 Anchor-2, LPP backbone)
%   OUT.overall  : summary struct (motif agreement + retention + magnitude deltas)
%   OUT.motifs   : extracted motif IDs (base vs sens)
%   OUT.info     : labels, alpha, sources, export paths

p = inputParser;
p.FunctionName = mfilename;

addParameter(p, 'BaseLabel', "baseline");
addParameter(p, 'SensLabel', "sensitivity");
addParameter(p, 'Alpha', 0.05);

addParameter(p, 'WriteCSV', true);
addParameter(p, 'WriteMAT', true);

addParameter(p, 'OutDir', "");
addParameter(p, 'OutStem', "");

addParameter(p, 'OutMAT', "");
addParameter(p, 'OutDetailCSV', "");
addParameter(p, 'OutOverallCSV', "");

parse(p, varargin{:});

baseLabel = string(p.Results.BaseLabel);
sensLabel = string(p.Results.SensLabel);
alpha = double(p.Results.Alpha);

writeCSV = logical(p.Results.WriteCSV);
writeMAT = logical(p.Results.WriteMAT);

Tb = as_table(base_csv_or_table);
Ts = as_table(sens_csv_or_table);

Tb = normalize_predictor_names(Tb);
Ts = normalize_predictor_names(Ts);

keys = ["N200 Anchor-2", "P300 Anchor-2", "LPP backbone"];

% Use a cell array so we can append structs safely even when first row differs
rows = {};

for i = 1: numel(keys)

    key = keys(i);

    rb = pick_row(Tb, key);
    rs = pick_row(Ts, key);

    R = struct();
    R.PredictorKey = string(key);

    % Motif IDs (if available)
    R.N200_Anchor2_MapID_Base = safe_col_value(rb, "N200_Anchor2_MapID");
    R.N200_Anchor2_MapID_Sens = safe_col_value(rs, "N200_Anchor2_MapID");

    R.P300_Anchor2_MapID_Base = safe_col_value(rb, "P300_Anchor2_MapID");
    R.P300_Anchor2_MapID_Sens = safe_col_value(rs, "P300_Anchor2_MapID");

    R.LPP_Backbone_Pair_Base = safe_col_value(rb, "LPP_Backbone_Pair");
    R.LPP_Backbone_Pair_Sens = safe_col_value(rs, "LPP_Backbone_Pair");

    R.LPP_Backbone_Pair_Base_Canon = canonical_pair(R.LPP_Backbone_Pair_Base);
    R.LPP_Backbone_Pair_Sens_Canon = canonical_pair(R.LPP_Backbone_Pair_Sens);

    % Compare motifs only if comparable (prevents inflated "same" counts when IDs are missing in that row)
    [R.MotifSame, R.MotifComparable] = motif_same_and_comparable(R.PredictorKey, R);

    % Interaction (confirmatory: within × valence)
    [R.Beta_Int_Base, R.p_Int_Base, R.Sig_Int_Base] = get_term(rb, "Beta", "pValue", alpha);
    [R.Beta_Int_Sens, R.p_Int_Sens, R.Sig_Int_Sens] = get_term(rs, "Beta", "pValue", alpha);
    [R.DirRet_Int, R.SigRet_Int] = retention(R.Beta_Int_Base, R.Sig_Int_Base, R.Beta_Int_Sens, R.Sig_Int_Sens);

    % Neutral slope
    [R.Beta_Neu_Base, R.p_Neu_Base, R.Sig_Neu_Base] = get_term(rb, "Beta_NeutralSlope", "pValue_NeutralSlope", alpha);
    [R.Beta_Neu_Sens, R.p_Neu_Sens, R.Sig_Neu_Sens] = get_term(rs, "Beta_NeutralSlope", "pValue_NeutralSlope", alpha);
    [R.DirRet_Neu, R.SigRet_Neu] = retention(R.Beta_Neu_Base, R.Sig_Neu_Base, R.Beta_Neu_Sens, R.Sig_Neu_Sens);

    % Valenced slope
    [R.Beta_Val_Base, R.p_Val_Base, R.Sig_Val_Base] = get_term(rb, "Beta_ValencedSlope", "pValue_ValencedSlope", alpha);
    [R.Beta_Val_Sens, R.p_Val_Sens, R.Sig_Val_Sens] = get_term(rs, "Beta_ValencedSlope", "pValue_ValencedSlope", alpha);
    [R.DirRet_Val, R.SigRet_Val] = retention(R.Beta_Val_Base, R.Sig_Val_Base, R.Beta_Val_Sens, R.Sig_Val_Sens);

    % Magnitude deltas
    R.AbsDelta_Int = abs(R.Beta_Int_Sens - R.Beta_Int_Base);
    R.AbsDelta_Neu = abs(R.Beta_Neu_Sens - R.Beta_Neu_Base);
    R.AbsDelta_Val = abs(R.Beta_Val_Sens - R.Beta_Val_Base);

    rows{end + 1, 1} = R;
end

if isempty(rows)
    detail = table();
else
    detail = struct2table([rows{:}]);
end

% Motif bundle (compact)
motifs = struct();
motifs.base = struct();
motifs.sens = struct();

motifs.base.Label = baseLabel;
motifs.sens.Label = sensLabel;

if isempty(detail)
    motifs.base.N200_Anchor2_MapID = NaN;
    motifs.sens.N200_Anchor2_MapID = NaN;
    motifs.base.P300_Anchor2_MapID = NaN;
    motifs.sens.P300_Anchor2_MapID = NaN;
    motifs.base.LPP_Backbone_Pair = "";
    motifs.sens.LPP_Backbone_Pair = "";
else
    motifs.base.N200_Anchor2_MapID = pick_one(detail.N200_Anchor2_MapID_Base);
    motifs.sens.N200_Anchor2_MapID = pick_one(detail.N200_Anchor2_MapID_Sens);

    motifs.base.P300_Anchor2_MapID = pick_one(detail.P300_Anchor2_MapID_Base);
    motifs.sens.P300_Anchor2_MapID = pick_one(detail.P300_Anchor2_MapID_Sens);

    motifs.base.LPP_Backbone_Pair = pick_one(detail.LPP_Backbone_Pair_Base_Canon);
    motifs.sens.LPP_Backbone_Pair = pick_one(detail.LPP_Backbone_Pair_Sens_Canon);
end

% Overall summary
overall = struct();
overall.Alpha = alpha;
overall.BaseLabel = baseLabel;
overall.SensLabel = sensLabel;

if isempty(detail)

    overall.MotifsCompared = 0;
    overall.MotifsSame_Count = 0;
    overall.MotifsAllSame = true;

    overall.BaselineSig_Int = 0;
    overall.DirRet_Int = 0;
    overall.SigRet_Int = 0;

    overall.BaselineSig_Neu = 0;
    overall.DirRet_Neu = 0;
    overall.SigRet_Neu = 0;

    overall.BaselineSig_Val = 0;
    overall.DirRet_Val = 0;
    overall.SigRet_Val = 0;

    overall.MedAbsDelta_Int = NaN;
    overall.MedAbsDelta_Neu = NaN;
    overall.MedAbsDelta_Val = NaN;

else

    mc = logical(detail.MotifComparable);
    ms = logical(detail.MotifSame);

    overall.MotifsCompared = nnz(mc);
    overall.MotifsSame_Count = nnz(ms & mc);

    if overall.MotifsCompared == 0
        overall.MotifsAllSame = true;
    else
        overall.MotifsAllSame = all(ms(mc));
    end

    overall.BaselineSig_Int = nnz(detail.Sig_Int_Base);
    overall.DirRet_Int = nnz(detail.DirRet_Int & detail.Sig_Int_Base);
    overall.SigRet_Int = nnz(detail.SigRet_Int & detail.Sig_Int_Base);

    overall.BaselineSig_Neu = nnz(detail.Sig_Neu_Base);
    overall.DirRet_Neu = nnz(detail.DirRet_Neu & detail.Sig_Neu_Base);
    overall.SigRet_Neu = nnz(detail.SigRet_Neu & detail.Sig_Neu_Base);

    overall.BaselineSig_Val = nnz(detail.Sig_Val_Base);
    overall.DirRet_Val = nnz(detail.DirRet_Val & detail.Sig_Val_Base);
    overall.SigRet_Val = nnz(detail.SigRet_Val & detail.Sig_Val_Base);

    overall.MedAbsDelta_Int = median(detail.AbsDelta_Int, 'omitnan');
    overall.MedAbsDelta_Neu = median(detail.AbsDelta_Neu, 'omitnan');
    overall.MedAbsDelta_Val = median(detail.AbsDelta_Val, 'omitnan');

end

OUT = struct();
OUT.detail = detail;
OUT.overall = overall;
OUT.motifs = motifs;

OUT.info = struct();
OUT.info.Alpha = alpha;
OUT.info.BaseLabel = baseLabel;
OUT.info.SensLabel = sensLabel;
OUT.info.BaseSource = source_tag(base_csv_or_table);
OUT.info.SensSource = source_tag(sens_csv_or_table);

% Decide export paths
outDir = string(p.Results.OutDir);
outStem = string(p.Results.OutStem);

if strlength(outDir) == 0
    outDir = infer_outdir(sens_csv_or_table);
end

if strlength(outStem) == 0
    outStem = "Q3_robustness_" + sanitize_label(baseLabel) + "_vs_" + sanitize_label(sensLabel);
end

outMAT = string(p.Results.OutMAT);
outDetailCSV = string(p.Results.OutDetailCSV);
outOverallCSV = string(p.Results.OutOverallCSV);

if strlength(outMAT) == 0
    outMAT = fullfile(outDir, outStem + ".mat");
end

if strlength(outDetailCSV) == 0
    outDetailCSV = fullfile(outDir, outStem + "_detail.csv");
end

if strlength(outOverallCSV) == 0
    outOverallCSV = fullfile(outDir, outStem + "_overall.csv");
end

OUT.info.OutDir = outDir;
OUT.info.OutStem = outStem;
OUT.info.OutMAT = outMAT;
OUT.info.OutDetailCSV = outDetailCSV;
OUT.info.OutOverallCSV = outOverallCSV;

% Export
if writeMAT
    try
        ensure_dir(outDir);
        save(outMAT, 'OUT');
    catch
    end
end

if writeCSV
    try
        ensure_dir(outDir);
        writetable(detail, outDetailCSV);

        Toverall = struct2table(overall, 'AsArray', true);
        writetable(Toverall, outOverallCSV);
    catch
    end
end

end


% ========================= HELPERS (local) =========================

function ensure_dir(d)

if ~exist(d, 'dir')
    mkdir(d);
end

end


function outDir = infer_outdir(x)

outDir = string(pwd);

if ischar(x) || isstring(x)
    fp = string(x);
    if strlength(fp) > 0
        try
            outDir = string(fileparts(fp));
            if strlength(outDir) == 0
                outDir = string(pwd);
            end
        catch
            outDir = string(pwd);
        end
    end
end

end


function s = sanitize_label(s)

s = string(s);
s = replace(s, " ", "_");
s = replace(s, "-", "_");
s = replace(s, ".", "_");
s = regexprep(s, '[^\w]', '_');
s = regexprep(s, '_+', '_');
s = strip(s, '_');

if strlength(s) == 0
    s = "label";
end

end


function tag = source_tag(x)

if istable(x)
    tag = "table_input";
elseif ischar(x) || isstring(x)
    tag = string(x);
else
    tag = "unknown_input";
end

end


function T = as_table(x)

if istable(x)
    T = x;
    return;
end

if ischar(x) || isstring(x)
    T = readtable(x);
    return;
end

error('Input must be a table or a CSV path.');
end


function T = normalize_predictor_names(T)

if ~ismember('Predictor', T.Properties.VariableNames)
    error('Expected a Predictor column in the Q3 output table.');
end

T.Predictor = string(T.Predictor);

T.Predictor = replace(T.Predictor, "Anchor 2", "Anchor-2");
T.Predictor = replace(T.Predictor, "Backbone", "backbone");

end


function row = pick_row(T, key)

row = table();
pat = string(key);

if pat == "N200 Anchor-2"
    mask = contains(T.Predictor, "N200", 'IgnoreCase', true) & contains(T.Predictor, "Anchor-2", 'IgnoreCase', true);
elseif pat == "P300 Anchor-2"
    mask = contains(T.Predictor, "P300", 'IgnoreCase', true) & contains(T.Predictor, "Anchor-2", 'IgnoreCase', true);
else
    mask = contains(T.Predictor, "LPP", 'IgnoreCase', true) & contains(T.Predictor, "backbone", 'IgnoreCase', true);
end

idx = find(mask, 1, 'first');

if isempty(idx)
    return;
end

row = T(idx, :);

end


function v = safe_col_value(Trow, colName)

v = NaN;

if isempty(Trow)
    return;
end

if ~ismember(colName, Trow.Properties.VariableNames)
    return;
end

x = Trow.(colName);

if iscell(x)
    x = x{1};
end

if isstring(x) || ischar(x)
    v = string(x);
elseif isnumeric(x)
    v = double(x(1));
else
    try
        v = string(x(1));
    catch
        v = NaN;
    end
end

end


function s = canonical_pair(x)

s = "";

if ~(isstring(x) || ischar(x))
    return;
end

str = strtrim(string(x));

if strlength(str) == 0
    return;
end

str = replace(str, "↔", "-");
str = replace(str, "<->", "-");
str = replace(str, "—", "-");
str = replace(str, "–", "-");
str = replace(str, ",", "-");

nums = regexp(str, '\d+', 'match');

if numel(nums) < 2
    return;
end

a = str2double(nums{1});
b = str2double(nums{2});

if ~isfinite(a) || ~isfinite(b)
    return;
end

ij = sort([a b]);
s = sprintf('%d-%d', ij(1), ij(2));

end


function [same, comparable] = motif_same_and_comparable(key, R)

same = false;
comparable = false;

if key == "N200 Anchor-2"

    a = force_scalar_double(R.N200_Anchor2_MapID_Base);
    b = force_scalar_double(R.N200_Anchor2_MapID_Sens);

    comparable = isfinite(a) && isfinite(b);
    same = comparable && (a == b);

elseif key == "P300 Anchor-2"

    a = force_scalar_double(R.P300_Anchor2_MapID_Base);
    b = force_scalar_double(R.P300_Anchor2_MapID_Sens);

    comparable = isfinite(a) && isfinite(b);
    same = comparable && (a == b);

else

    b = string(R.LPP_Backbone_Pair_Base_Canon);
    s = string(R.LPP_Backbone_Pair_Sens_Canon);

    comparable = (strlength(b) > 0) && (strlength(s) > 0);
    same = comparable && (b == s);

end

end


function x = force_scalar_double(v)

x = NaN;

try
    if isnumeric(v)
        x = double(v(1));
        return;
    end

    if isstring(v) || ischar(v)
        x = str2double(string(v));
        return;
    end

    if iscell(v) && ~isempty(v)
        x = str2double(string(v{1}));
        return;
    end
catch
    x = NaN;
end

end


function [beta, p, sig] = get_term(Trow, betaCol, pCol, alpha)

beta = NaN;
p = NaN;
sig = false;

if isempty(Trow)
    return;
end

if ismember(betaCol, Trow.Properties.VariableNames)
    beta = force_scalar_double(Trow.(betaCol));
end

if ismember(pCol, Trow.Properties.VariableNames)
    p = force_scalar_double(Trow.(pCol));
end

sig = isfinite(p) && p < alpha;

end


function [dirRet, sigRet] = retention(betaB, sigB, betaS, sigS)

dirRet = false;
sigRet = false;

if ~sigB
    return;
end

if isfinite(betaB) && isfinite(betaS)
    dirRet = sign(betaB) == sign(betaS);
end

sigRet = sigS;

end


function v = pick_one(x)

v = NaN;

if isempty(x)
    return;
end

try
    if iscell(x)
        x = string(x);
    end

    if isstring(x)
        idx = find(strlength(x) > 0 & x ~= "NaN", 1, 'first');
        if ~isempty(idx)
            v = x(idx);
        end
    else
        idx = find(isfinite(x), 1, 'first');
        if ~isempty(idx)
            v = x(idx);
        end
    end
catch
    v = NaN;
end

end
