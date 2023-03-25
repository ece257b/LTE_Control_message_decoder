function [conflict_flag, new_potential_rnti, new_conflict_rnti] = pdcch_filter(enb, rxsubframe, potential_rnti, hest, nest, conflict_flag, conflict_rnti)
    conflict_flag = false;
    new_potential_rnti = dictionary;
    new_conflict_rnti = {};
    for i = 1:numel(conflict_rnti)
        if(isKey(potential_rnti,conflict_rnti{i}.left) && isKey(potential_rnti, conflict_rnti{i}.right))
            
            pdsch = setPDSCHConfiguration(enb, potential_rnti(conflict_rnti{i}.left).dcimsg, conflict_rnti{i}.left);
            if ~isempty(pdsch)
                pdsch.NTurboDecIts = 5;
                pdschIndices = ltePDSCHIndices(enb, pdsch, pdsch.PRBSet);
                [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsubframe, hest);
                dlschBits = ltePDSCHDecode(enb, pdsch, pdschRx, pdschHest, nest);
                if(all(abs(dlschBits{1,1}) > 1e-6))
                    new_potential_rnti(conflict_rnti{i}.left) = potential_rnti(conflict_rnti{i}.left);
                end
            else
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.right);
            end
            
            tmp=potential_rnti(conflict_rnti{i}.right);
            pdsch = setPDSCHConfiguration(enb, potential_rnti(conflict_rnti{i}.right).dcimsg, conflict_rnti{i}.right);
            if ~isempty(pdsch)
                pdsch.NTurboDecIts = 5;
                pdschIndices = ltePDSCHIndices(enb, pdsch, pdsch.PRBSet);
                [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsubframe, hest);
                dlschBits = ltePDSCHDecode(enb, pdsch, pdschRx, pdschHest, nest);
                if(all(abs(dlschBits{1,1}) > 1e-6))
                    new_potential_rnti(conflict_rnti{i}.right) = potential_rnti(conflict_rnti{i}.right);
                end
            else
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.right);
            end
        end
    end
    if (numEntries(new_potential_rnti) > 1)
        potential_rnti_list = keys(new_potential_rnti);
        for m = 1:length(potential_rnti_list)
           for n = 2:length(potential_rnti_list)
               vecLoc1 = potential_rnti(potential_rnti_list(m)).Loc;
               vecLoc2 = potential_rnti(potential_rnti_list(n)).Loc;
               if(sum(ismember(vecLoc1,vecLoc2)) > 0)
                   conflict_rntis = struct;
                   conflict_rntis.left = new_potential_rnti_list(m);
                   conflict_rntis.right = new_potential_rnti_list(n);
                   new_conflict_rnti{idx} = conflict_rntis;
                   idx = idx+1;
                   conflict_flag = true;
               end
           end
        end
    end
end