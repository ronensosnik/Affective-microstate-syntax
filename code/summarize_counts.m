function summary = summarize_counts(detail, ERPs)
% Build a compact summary table: overall and by ERP

overall = table;
overall.Scope = "Overall";
overall.Total_Models        = height(detail);
overall.Any_FDR_Change      = sum(detail.Any_FDR_Change);
overall.Direction_Changed   = sum(detail.Direction_Changed);

overall.Lost_Any = sum(detail.Lost_Int | detail.Lost_Grp | detail.Lost_Cond);
overall.Gained_Any = sum(detail.Gained_Int | detail.Gained_Grp | detail.Gained_Cond);

% Per-ERP rows
rows = overall;
for e = 1:numel(ERPs)
    sub = detail(strcmp(detail.ERP, ERPs{e}), :);
    t = table;
    t.Scope = string(ERPs{e});
    t.Total_Models      = height(sub);
    t.Any_FDR_Change    = sum(sub.Any_FDR_Change);
    t.Direction_Changed = sum(sub.Direction_Changed);
    t.Lost_Any          = sum(sub.Lost_Int | sub.Lost_Grp | sub.Lost_Cond);
    t.Gained_Any        = sum(sub.Gained_Int | sub.Gained_Grp | sub.Gained_Cond);
    rows = [rows; t]; %#ok<AGROW>
end

summary = rows;
end