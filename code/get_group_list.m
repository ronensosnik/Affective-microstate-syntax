function list = get_group_list(France, grp)
switch grp
    case 'BP_I_Depressed'
        list = France.BP_I.Depressed;
    case 'BP_I_Euthymic'
        list = France.BP_I.Euthymic;
    case 'BP_II_Depressed'
        list = France.BP_II.Depressed;
    case 'BP_II_Euthymic'
        list = France.BP_II.Euthymic;
    case 'HC'
        list = France.HC;
    case 'Siblings'
        list = France.Siblings;
    otherwise
        error('Unknown group: %s', grp);
end
end
