function [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx)
P = (S.C_segments(:) / S.M_segments);     % 7×1
E = S.N_transitions * (P * P.');          % 7×7 independence
yvec = S.trans_counts(sub2ind([7,7], srcIdx, tgtIdx));
evec = E(sub2ind([7,7], srcIdx, tgtIdx));
end