function [conflict_flag, potential_rnti] = Active_freq_filter(conflict_rnti, Active_Rnti, potential_rnti)
    for i = 1:numel(conflict_rnti)
        if(isKey(potential_rnti,conflict_rnti{i}.left) && isKey(potential_rnti, conflict_rnti{i}.right))
            lefttime = get(Active_Rnti,num2str(conflict_rnti{i}.left));
            righttime = get(Active_Rnti,num2str(conflict_rnti{i}.right));
            if(lefttime >= righttime)
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.right);
            else
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.left);
            end
        end
    end
    conflict_flag = false;
end