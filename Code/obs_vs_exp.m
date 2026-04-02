function [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx)

%  The below code does not take in consideration and omits self transitions (e.g., A -> A) so the expected count is smaller and it can lead to more  transitions than deviate from independence

if 0
    P = (S.C_segments(:) / S.M_segments); % 7×1
    E = S.N_transitions * (P * P.'); % 7×7 independence
    yvec = S.trans_counts(sub2ind([7, 7], srcIdx, tgtIdx));
    evec = E(sub2ind([7, 7], srcIdx, tgtIdx));
    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
end

%  The below code DOES take in consideration and omits self transitions (e.g., A -> A)

if 1
    P = (S.C_segments(:) / S.M_segments); % 7×1
    Ssq = sum(P.^2);
    den = 1 - Ssq;
    kK=numel(S.C_segments);

    gamma = 1 / den;

    E = gamma * S.N_transitions * (P * P.');
    E(1: kK + 1: kK * kK) = NaN;

    yvec = S.trans_counts(sub2ind([kK, kK], srcIdx, tgtIdx));
    evec = E(sub2ind([kK, kK], srcIdx, tgtIdx));

    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
end

end