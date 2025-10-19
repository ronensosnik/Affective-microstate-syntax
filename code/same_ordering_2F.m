function tf = same_ordering_2F(Tb, Tm, withinF, byF)
% Compare ordering of EMM_BT across levels of withinF, within each level of byF.
% Returns true if orderings match for all slices.
tf = true;
if isempty(Tb) || isempty(Tm), tf = true; return; end

lvW = categories(Tb.(withinF));
lvB = categories(Tb.(byF));

% Ensure Tm has same factor levels (otherwise, assume mismatch)
if ~all(ismember(lvW, categories(Tm.(withinF)))) || ~all(ismember(lvB, categories(Tm.(byF))))
    tf = false; return;
end

for j = 1:numel(lvB)
    b = lvB{j};
    xb = Tb.EMM_BT(Tb.(byF)==b);
    gw = Tb.(withinF)(Tb.(byF)==b);

    xm = Tm.EMM_BT(Tm.(byF)==b);
    gm = Tm.(withinF)(Tm.(byF)==b);

    % Align by level names
    [gw_u, ib] = sortrows(cellstr(gw));
    [gm_u, im] = sortrows(cellstr(gm));

    if ~isequal(gw_u, gm_u)
        tf = false; return;
    end

    xb = xb(ib);
    xm = xm(im);

    % Compare rank order via sign of successive differences
    sb = sign(diff(xb));
    sm = sign(diff(xm));
    if numel(sb) ~= numel(sm) || any(sb ~= sm)
        tf = false; return;
    end
end
end
