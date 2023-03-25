function [conflict_flag, potential_rnti, new_conflict_rnti] = format_filter(enbConfig,ueConfig, potential_rnti, conflict_rnti)
    enbConfig = mwltelibrary('validateLTEParameters',enbConfig,'CellRefP','mandatory');
    ueConfig  = mwltelibrary('validateLTEParameters',ueConfig, 'EnableCarrierIndication',...
                                                               'EnableSRSRequest',...
                                                               'EnableMultipleCSIRequest',...
                                                               'NTxAnts');
    conflict_flag = false;
    pdcchConfig = ueConfig;
    new_conflict_rnti = {};
    pdcchConfig.ControlChannelType = 'PDCCH';
    idx = 1;
    if isfield(pdcchConfig,'DCIFormat')
        pdcchConfig = rmfield(pdcchConfig,'DCIFormat');
    end 
    pdcchConfig.SearchSpace = 'UESpecific';
    if(pdcchConfig.SearchSpace == "UESpecific")         
        startingPdcchFormat = 0;
    else % common search space
        startingPdcchFormat = 2;
    end
    for i = 1:numel(conflict_rnti)
        resolveflag = 0;
        if(isKey(potential_rnti,conflict_rnti{i}.left))
            pdcchFormat = log2(length(potential_rnti(conflict_rnti{i}.left).Loc));
            pdcchConfig.PDCCHFormat = pdcchFormat;
            pdcchConfig.RNTI=conflict_rnti{i}.left;
            pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});
            pdcchCandidates = unique(pdcchCandidates,'rows');
            pdcchCandidatesDims = size(pdcchCandidates);
            noOfCandidates = pdcchCandidatesDims(1);
            Locflag = false;
            for candidate=1:noOfCandidates
                tmp = potential_rnti(conflict_rnti{i}.left).Loc;
                if(isequal(potential_rnti(conflict_rnti{i}.left).Loc, (pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,2)/72)-1)))
                    resolveflag = resolveflag+1;
                    Locflag = true;
                    break;
                end
            end
            if(Locflag == false)
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.left);
            end
        end
        if(isKey(potential_rnti,conflict_rnti{i}.right))
            pdcchFormat = log2(length(potential_rnti(conflict_rnti{i}.right).Loc));
            pdcchConfig.PDCCHFormat = pdcchFormat;
            pdcchConfig.RNTI=conflict_rnti{i}.left;
            pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});
            pdcchCandidates = unique(pdcchCandidates,'rows');
            pdcchCandidatesDims = size(pdcchCandidates);
            noOfCandidates = pdcchCandidatesDims(1);
            Locflag = false;
            for candidate=1:noOfCandidates
                if(isequal(potential_rnti(conflict_rnti{i}.right).Loc, (pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,2)/72)-1)))
                    resolveflag = resolveflag+1;
                    Locflag = true;
                    break;
                end
            end
            if(Locflag == false)
                potential_rnti=delete_key(potential_rnti,conflict_rnti{i}.right);
            end
        end
        if(resolveflag == 2)
            new_conflict_rnti{idx} = conflict_rnti{i};
            conflict_flag = true;
            idx = idx+1;
        end
    end
end