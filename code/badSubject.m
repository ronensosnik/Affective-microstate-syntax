function tf = badSubject(S)
tf = isempty(S) || ~isfield(S,'trans_counts') || isempty(S.trans_counts) || S.M_segments<=0 || S.N_transitions<=0;
end