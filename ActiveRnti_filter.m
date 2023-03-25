function [conflict_flag, potential_rnti, conflict_rnti] = ActiveRnti_filter(rnti_dict, Active_Rnti)
    conflict_flag = false;
    potential_rnti = dictionary;
    conflict_rnti = {};
    idx = 1;
    if (numEntries(rnti_dict.rawrnti) > 0)
        rnti_list = keys(rnti_dict.rawrnti);
        for j = 1:length(rnti_list)
            if(iscached(Active_Rnti,num2str(rnti_list(j))))
                potential_rnti(rnti_list(j)) = rnti_dict.rawrnti(rnti_list(j));
            end
        end
        if (numEntries(potential_rnti) > 1)
            potential_rnti_list = keys(potential_rnti);
            for m = 1:length(potential_rnti_list)
               for n = 2:length(potential_rnti_list)
                   vecLoc1 = potential_rnti(potential_rnti_list(m)).Loc;
                   vecLoc2 = potential_rnti(potential_rnti_list(n)).Loc;
                   if(sum(ismember(vecLoc1,vecLoc2)) > 0)
                       conflict_rntis = struct;
                       conflict_rntis.left = potential_rnti_list(m);
                       conflict_rntis.right = potential_rnti_list(n);
                       conflict_rnti{idx} = conflict_rntis;
                       idx = idx+1;
                       conflict_flag = true;
                   end
               end
            end
        end
    end
end