function MedMap = build_medmap(France, Groups)
MedMap = struct();
for g = 1:numel(Groups)
    grp = Groups{g};
    list = get_group_list(France, grp);
    n = numel(list);
    AD  = zeros(n,1); AP  = zeros(n,1); MS = zeros(n,1);
    ANX = zeros(n,1); OTH = zeros(n,1); PC = zeros(n,1);
    for i = 1:n
        s = list{i};
        AD(i)  = getfield_safe(s,'AD',0);
        AP(i)  = getfield_safe(s,'AP',0);
        MS(i)  = getfield_safe(s,'MS',0);
        ANX(i) = getfield_safe(s,'ANX',0);
        OTH(i) = getfield_safe(s,'OTHER',0);
        PC(i)  = getfield_safe(s,'poly_count',0);
    end
    MedMap.(grp).AD = AD; MedMap.(grp).AP = AP; MedMap.(grp).MS = MS;
    MedMap.(grp).ANX = ANX; MedMap.(grp).OTHER = OTH; MedMap.(grp).poly_count = PC;
end
end