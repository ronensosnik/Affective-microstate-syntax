function C = transitionCounts(MSClass, win)
% 7×7 directed transition counts within samples win(1):win(2)

C = zeros(7, 7);

t1 = win(1);
t2 = win(2);

if t1 < 1
    t1 = 1;
end

if t2 > size(MSClass, 1)
    t2 = size(MSClass, 1);
end

for tr = 1: size(MSClass, 2)
    l = MSClass(t1: t2, tr);
    l = double(l);
    valid = isfinite(l) & l >= 1 & l <= 7;
    l(~valid) = 0;

    pre = l(1: end - 1);
    nxt = l(2: end);
    mask = pre ~= nxt & pre > 0 & nxt > 0;

    pre = pre(mask);
    nxt = nxt(mask);

    if isempty(pre)
        continue
    end

    for k = 1: numel(pre)
        C(pre(k), nxt(k)) = C(pre(k), nxt(k)) + 1;
    end
end

end