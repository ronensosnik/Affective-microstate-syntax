function C = readme_lines()

C = {
    'Metrics_Meds_Robustness.xlsx';
    'Sheets follow the pattern: <Measure>_<Variant>_<detail|summary>';
    '  Measure: Coverage | Duration | Occurrence';
    '  Variant: classes (AD,AP,MS,ANX,OTHER) | poly (poly_count only)';
    'detail: 21 rows (3 ERPs x 7 Maps), flags for FDR changes and direction flips';
    'summary: overall + per-ERP counts (Any_FDR_Change, Lost_Any, Gained_Any, Direction_Changed)';
    'All categoricals converted to strings for portability.';
    };
end
