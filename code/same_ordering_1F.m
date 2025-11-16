function tf = same_ordering_1F(Tb, Tm, fac)

% Compare ordering of EMM_BT across levels of one factor

tf = true;

if isempty(Tb) || isempty(Tm)
    tf = true;
    return;
end

lv = categories(Tb.(fac));

if ~all(ismember(lv, categories(Tm.(fac))))
    tf = false;
    return;
end

% Align by level names

[lvb, ib] = sortrows(cellstr(Tb.(fac)));
[lvm, im] = sortrows(cellstr(Tm.(fac)));

if ~isequal(lvb, lvm)
    tf = false;
    return;
end

xb = Tb.EMM_BT(ib);
xm = Tm.EMM_BT(im);

sb = sign(diff(xb));
sm = sign(diff(xm));

if numel(sb) ~= numel(sm) || any(sb ~= sm)
    tf = false;
end
end