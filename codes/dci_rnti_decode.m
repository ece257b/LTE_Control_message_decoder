
function [decRnti_dict_mat, candidate_Rnti, Active_Rnti] = dci_rnti_decode(enbConfig,ueConfig,inBits, candidate_Rnti, Active_Rnti)

    % Validate and default the shared (mandatory and optional) parameters
    enbConfig = mwltelibrary('validateLTEParameters',enbConfig,'CellRefP','mandatory');
    ueConfig  = mwltelibrary('validateLTEParameters',ueConfig, 'EnableCarrierIndication',...
                                                               'EnableSRSRequest',...
                                                               'EnableMultipleCSIRequest',...
                                                               'NTxAnts');
    % Input size validation
    if length(inBits)<72
        error('lte:error','Input soft bits should be at least a CCE in length (72 bits).');
    end
    
    % Deduce total number of REGs (1 REG = 4 REs = 8 bits)
    enbConfig.NREG = floor(length(inBits)/(72/9));
    
    % DCI formats for common and UE-specific search space
    %
    % Common search space DCI formats, note that Format1A is not listed
    % because it is the same size as Format0
    dciFormats{1} = {'Format0','Format1C'};
  
    % Cache the UE specific parameters and initially remove the RNTI since 
    % we will start by searching in the common search space
    pdcchConfig = ueConfig;
    pdcchConfig.ControlChannelType = 'PDCCH';
    if isfield(pdcchConfig,'RNTI')
        pdcchConfig = rmfield(pdcchConfig,'RNTI');
    end
    % Ensure that 'DCIFormat' is not present in UE-specific parameters
    % since it may conflict with the use of lteDCI later 
    if isfield(pdcchConfig,'DCIFormat')
        pdcchConfig = rmfield(pdcchConfig,'DCIFormat');
    end   
    % Get DCI formats and lengths (relative to common search space 
    % for formats 0/1A and UE-specific otherwise), only the first format
    % for each unique message length is listed
    pdcchConfig.SearchSpace = 'UESpecific';
    dciinfo = lte.internal.uniqueDCILengths(enbConfig,pdcchConfig);
     
    % Identify UE-specific search space DCI formats
    dcimessages = fieldnames(dciinfo);
    uespecific_exclusionList = {'Format3','Format3A','Format0','Format4'};   % Remove 3/3A from this list
    dciFormats{2} = setdiff(dcimessages,uespecific_exclusionList);
       
    % Intermediate local variables
    idx = 1;
    decDCI = {};
    decDCIBits = {};
    reservedLoc = [];
    possible_rnti = dictionary;
    decRnti_dict_mat = dictionary;

    if(pdcchConfig.SearchSpace == "UESpecific")         
        % PDCCH format for UE-specific search space can be 0,1,2 or 3
        startingPdcchFormat = 0;
        searchType = 2;
    else % common search space
        % PDCCH format for common search space can either be 2 or 3 (i.e.
        % aggregation level of 4 or 8 CCEs)
        startingPdcchFormat = 2;
        searchType = 1;
    end
    % Get the message lengths for the current search space, only the 
    % first format for each unique message length is listed
    dciinfo = lte.internal.uniqueDCILengths(enbConfig,pdcchConfig);
    for pdcchFormat = 3:-1:startingPdcchFormat
        pdcchConfig.PDCCHFormat = pdcchFormat;
%         pdcchConfig.RNTI=255;

        % Performs common and/or ue-specific search depending on whether
        % the RNTI field is present in the channel/UE-specific parameters
        % TODO optimize by finding peak energy.
        pdcchCandidates = ltePDCCHSpace(enbConfig,pdcchConfig,{'bits','1based'});
        pdcchCandidates = generate_candidate(abs(pdcchCandidates(1,2)-pdcchCandidates(1,1)+1),length(inBits));
        % PDCCH candidates need not to be unique so picking the unique set
        % of candidates to optimize the search
%         pdcchCandidates = unique(pdcchCandidates,'rows');

        pdcchCandidatesDims = size(pdcchCandidates);
        noOfCandidates = pdcchCandidatesDims(1);
        
        % Test each candidate in the search space for the presence of a decodable message format
        % TODO: Deal with multiple dci in one PDCCH.
        for candidate=1:noOfCandidates
%             if(sum(reservedLoc == pdcchCandidates(candidate,1)/72) == 0)
            if((pdcchCandidates(candidate,1)<length(inBits)) && (pdcchCandidates(candidate,2)<=length(inBits)))
                input = inBits(pdcchCandidates(candidate,1):pdcchCandidates(candidate,2));
                CCE = reshape(input,72,[]);
                dcienergy = all(abs(CCE) > 1e-6);
                if all(dcienergy(:)>0)
                    for dciFormatIdx=1:length(dciFormats{searchType}) % Iterating through all DCI formats                     
                        df = dciFormats{searchType}{dciFormatIdx};        
                        mlength = dciinfo.(df);
                        [dciMessageBits,decRnti] = lteDCIDecode(mlength,input);
                        pdcchConfig.RNTI = decRnti;
                        [dciMessage, dciMessageBits] = lteDCI(enbConfig,pdcchConfig,dciMessageBits);
                        legal_rnti=Loc_filter(enbConfig,pdcchConfig,(pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,1)/72)+(2^pdcchFormat)-1));
                        if(legal_rnti)
                            rnti_dciinfo = struct;
                            rnti_dciinfo.dcimsg = dciMessage;
                            rnti_dciinfo.Loc = (pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,1)/72)+(2^pdcchFormat)-1);
                            if(numEntries(decRnti_dict_mat) > 0 && isKey(decRnti_dict_mat, decRnti) && sum(ismember(rnti_dciinfo.Loc, decRnti_dict_mat(decRnti).Loc))>0)
                                continue;
                            end
                            decRnti_dict_mat(decRnti) = rnti_dciinfo;
                            if candidate_Rnti.iscached(num2str(decRnti))
                                if Active_Rnti.iscached(num2str(decRnti))
                                    activetime = Active_Rnti.get(num2str(decRnti));
                                    Active_Rnti.put(num2str(decRnti),activetime+1);
                                else
                                    Active_Rnti.put(num2str(decRnti),1);
                                end
                            end
                            candidate_Rnti.put(num2str(decRnti),dciMessage.Allocation);
                        end
%                             if(ueConfig.RNTI == decRnti && ~any(reservedLoc == pdcchCandidates(candidate,1)/72))
%                                 % Creating DCI message for successfully decoded PDCCH payload bits
%                                 [dciMessage,dciMessageBits] = lteDCI(enbConfig,pdcchConfig,dciMessageBits);
%                                 decDCIBits{idx} = dciMessageBits;
%                                 decDCI{idx} = dciMessage;
%                                 reservedLoc = [reservedLoc,(pdcchCandidates(candidate,1)/72:(pdcchCandidates(candidate,1)/72)+(2^pdcchFormat)-1)]; %#ok<AGROW>
%                                 idx = idx + 1;
%                                 break;     % Now break out of the format list loop since a format was found in the candidate
%                             end
                    end
%                     end
                end
            end
        end
    end
end
