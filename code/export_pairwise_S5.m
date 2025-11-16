function export_pairwise_S5(OUT, measure, pipeline_label, out_csv)

% EXPORT_PAIRWISE_S5
% Flatten all BH–FDR-significant pairwise contrasts from your OUT struct
% into a single CSV (S5-style).
%
% Usage examples:
%   export_pairwise_S5(Out_Coverage_base,      'Coverage',   'baseline',   'S5_pairs_baseline_Coverage.csv');
%   export_pairwise_S5(Out_Coverage_sens_peaks,'Coverage',   'peaks_only', 'S5_pairs_peaks_Coverage.csv');
%   export_pairwise_S5(Out_Duration_base,      'Duration',   'baseline',   'S5_pairs_baseline_Duration.csv');
%   export_pairwise_S5(Out_Duration_sens_peaks,'Duration',   'peaks_only', 'S5_pairs_peaks_Duration.csv');
%   export_pairwise_S5(Out_Occurrence_base,    'Occurrence', 'baseline',   'S5_pairs_baseline_Occurrence.csv');
%   export_pairwise_S5(Out_Occurrence_sens_peaks,'Occurrence','peaks_only','S5_pairs_peaks_Occurrence.csv');

ERPs = {'N200', 'P300', 'LPP'};
Micros = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};

rows = {}; % will collect {Measure, ERP, Micro, Scope, Factor, Context, L1, L2, Eff, CIlo, CIhi, p, q, Units, Pipeline}

for e = 1: 3
    for m = 1: 7

        R = OUT.(measure){e,m};

        if isempty(R) || ~isfield(R,'pairwise') || isempty(R.pairwise)
            continue;
        end

        % ---- No-interaction main-effect pairwise blocks (most common) ----

        if isfield(R.pairwise,'Group') && ~isempty(R.pairwise.Group)
            rows = [rows; dumpPairTable(R.pairwise.Group,     measure, ERPs{e}, Micros{m}, 'Group (marginal)', 'Group', '', pipeline_label)];
        end

        if isfield(R.pairwise,'Condition') && ~isempty(R.pairwise.Condition)
            rows = [rows; dumpPairTable(R.pairwise.Condition, measure, ERPs{e}, Micros{m}, 'Condition (marginal)', 'Condition', '', pipeline_label)];
        end

        % ---- Simple-effects (only if you had interactions; harmless otherwise) ----

        if isfield(R.pairwise,'Group_within_Condition') && ~isempty(R.pairwise.Group_within_Condition)
            rows = [rows; dumpPairTable(R.pairwise.Group_within_Condition, measure, ERPs{e}, Micros{m}, 'Group @ Condition', 'Group', 'Condition', pipeline_label)];
        end

        if isfield(R.pairwise,'Condition_within_Group') && ~isempty(R.pairwise.Condition_within_Group)
            rows = [rows; dumpPairTable(R.pairwise.Condition_within_Group, measure, ERPs{e}, Micros{m}, 'Condition @ Group', 'Condition', 'Group', pipeline_label)];
        end
    end
end

if isempty(rows)
    warning('No BH–FDR-significant pairwise contrasts found.');
    T = cell2table(cell(0, 15), 'VariableNames', {'Measure','ERP','Microstate','Scope','Factor','Context','Level1','Level2','Effect','CI_low','CI_high','p','pFDR','Units','Pipeline'});
    writetable(T, out_csv);
    return;
end

T = cell2table(rows, 'VariableNames', {'Measure','ERP','Microstate','Scope','Factor','Context','Level1','Level2','Effect','CI_low','CI_high','p','pFDR','Units','Pipeline'});

writetable(T, out_csv);
fprintf('Wrote %d rows to %s\n', height(T), out_csv);
end