function legal_rnti = Loc_filter(enbConfig,ueConfig, Loc)
    enbConfig = mwltelibrary('validateLTEParameters',enbConfig,'CellRefP','mandatory');
    ueConfig  = mwltelibrary('validateLTEParameters',ueConfig, 'EnableCarrierIndication',...
                                                               'EnableSRSRequest',...
                                                               'EnableMultipleCSIRequest',...
                                                               'NTxAnts');
    pdcchConfig = ueConfig;
    pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});
    pdcchCandidates = unique(pdcchCandidates,'rows');
    pdcchCandidatesDims = size(pdcchCandidates);
    noOfCandidates = pdcchCandidatesDims(1);
    legal_rnti = false;
    for candidate=1:noOfCandidates
        if(isequal(Loc, (pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,2)/72)-1)))
            legal_rnti = true;
            break;
        end
    end
end