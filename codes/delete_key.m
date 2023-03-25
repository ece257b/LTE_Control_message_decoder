function new_potential_rnti=delete_key(potential_rnti,dictkey)
    new_potential_rnti = dictionary;
    rnti_list = keys(potential_rnti);
    for i = 1:length(rnti_list)
        if rnti_list(i) == dictkey
            continue;
        end
        new_potential_rnti(rnti_list(i)) = potential_rnti(rnti_list(i));
    end
end