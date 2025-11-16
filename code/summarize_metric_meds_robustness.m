function OUT = summarize_metric_meds_robustness(Out_base, Out_med, measure, alphaFDR)

% SUMMARIZE_METRIC_MEDS_ROBUSTNESS
% Compare baseline vs medication-adjusted LME results for a static metric.
%
% Usage:
%   OUT = summarize_metric_meds_robustness(Out_base, Out_med, 'Coverage', 0.05);
%
% Inputs:
%   Out_base : struct from lme_posthoc_21FDR_transformsBT (no meds)
%   Out_med  : struct from lme_posthoc_21FDR_transformsBT (with meds)
%   measure  : 'Duration' | 'Coverage' | 'Occurrence'
%   alphaFDR : FDR threshold (default 0.05)
%
% Output (struct OUT):
%   OUT.detail  : table, 21 rows (ERP×Map) with flags for FDR changes and direction flips
%   OUT.summary : table with overall counts and per-ERP counts

if nargin < 4 || isempty(alphaFDR)
 alphaFDR = 0.05; 
end

measure = validatestring(measure, {'Duration', 'Coverage', 'Occurrence'});

ERPs = Out_med.info.ERPs;
Micros = Out_med.info.Microstates;

rows = []; 
r = 0;

for e = 1: numel(ERPs)
    for m = 1: numel(Micros)

        r = r + 1;
        Rb = Out_base.(measure){e, m};
        Rm = Out_med.(measure){e, m};

        % --- FDR flags (baseline vs med-adjusted) ---

        bI = lt_alpha(safeget(Rb, 'tests', 'interaction', 'pFDR'), alphaFDR);
        bG = lt_alpha(safeget(Rb, 'tests', 'group', 'pFDR'), alphaFDR);
        bC = lt_alpha(safeget(Rb,'tests', 'condition', 'pFDR'), alphaFDR);

        mI = lt_alpha(safeget(Rm, 'tests', 'interaction', 'pFDR'), alphaFDR);
        mG = lt_alpha(safeget(Rm, 'tests', 'group', 'pFDR'), alphaFDR);
        mC = lt_alpha(safeget(Rm, 'tests', 'condition', 'pFDR'), alphaFDR);

        lostI = bI & ~mI; 
        gainedI = ~bI & mI;
        lostG = bG & ~mG; 
        gainedG = ~bG & mG;
        lostC = bC & ~mC; 
        gainedC = ~bC & mC;

        anyFDRchange = (lostI|gainedI) | (lostG|gainedG) | (lostC|gainedC);

        % --- Direction change (only where BOTH models are FDR-sig in the same test type) ---

        dirChanged = false;

        try
            if bI && mI

                % Compare Group-by-Condition EMM orderings

                EB = safeget(Rb, 'emm', 'Group_by_Condition');
                EM = safeget(Rm, 'emm', 'Group_by_Condition');
                dirChanged = ~same_ordering_2F(EB, EM, "Group", "Condition");

            elseif ~bI && ~mI && bG && mG
                EB = safeget(Rb, 'emm', 'Group');
                EM = safeget(Rm, 'emm', 'Group');
                dirChanged = ~same_ordering_1F(EB, EM, "Group");

            elseif ~bI && ~mI && bC && mC
                EB = safeget(Rb, 'emm', 'Condition');
                EM = safeget(Rm, 'emm', 'Condition');
                dirChanged = ~same_ordering_1F(EB, EM, "Condition");
                
            else
                dirChanged = false;
            end
        catch

            % If anything is missing, treat as no direction change

            dirChanged = false;
        end

        rows(r).ERP = ERPs{e};
        rows(r).Map = Micros{m};
        rows(r).Measure = measure;

        rows(r).Base_Int = bI; 
        rows(r).Med_Int = mI;  
        rows(r).Lost_Int = lostI; 
        rows(r).Gained_Int = gainedI;
        rows(r).Base_Grp = bG;
        rows(r).Med_Grp = mG;  
        rows(r).Lost_Grp = lostG;   
        rows(r).Gained_Grp = gainedG;
        rows(r).Base_Cond = bC;
        rows(r).Med_Cond = mC;
        rows(r).Lost_Cond = lostC;   
        rows(r).Gained_Cond = gainedC;

        rows(r).Any_FDR_Change = anyFDRchange;
        rows(r).Direction_Changed = dirChanged;
    end
end

detail = struct2table(rows);

% ---- Summary (overall + by ERP) ----

summary = summarize_counts(detail, ERPs);
OUT = struct('detail', detail, 'summary', summary);

end % ======= END MAIN =======