function rnti_decoder(varargin)
% clear;
% close all;

%% Param input check
if (nargin==0)
    filename='../lte_test_R2_2frame.dat';
    sample_rate = 15.36e6;
elseif (nargin==1)
    filename=varargin{1};
    sample_rate = 15.36e6;
elseif (nargin==2)
    filename = varargin{1};
    sample_rate = str2double(varargin{2});
else
    fprintf('\nOnly need filename and sample_rate');
    return;
end

candidate_Rnti = LRU(100);
Active_Rnti = LRU(10);

%% Load iq samples from file
fid = fopen(filename,'r');
data = fread(fid,[2,Inf], 'float');
Data = data(1,:) + data(2,:)*1i;
[r,c] = size(Data);
Data = reshape(Data, c, r);
if (isempty(Data))
    fprintf('\nReceived signal is empty.\n');
    return;
end


%% Downsample data for cell search
% Decode primary and secondary synchronization signals (PSS and SSS) as well as PBCH (contains MIB).
% The PBCH and PSS/SSS information all lie in the central 6 resource
% blocks. Resample iq waveform accordingly. 

enb = struct;
enb.NDLRB = 6;
enb.CyclicPrefix = 'Normal';
ofdmInfo = lteOFDMInfo(enb); %Set up ofdm information for demodulation

if (sample_rate~=ofdmInfo.SamplingRate)
    if (sample_rate < ofdmInfo.SamplingRate)
        warning('The received signal sampling rate (%0.3fMs/s) is lower than the desired sampling rate for cell search / RNTI decoding (%0.3fMs/s); cell search / RNTI decoding may fail.',sameple_rate/1e6,ofdmInfo.SamplingRate/1e6);
    end
    fprintf('\nResampling from %0.3fMs/s to %0.3fMs/s for cell search / RNTI decoding...\n',sample_rate/1e6,ofdmInfo.SamplingRate/1e6);
else
    fprintf('\nResampling not required; received signal is at desired sampling rate for cell search / RNTI decoding (%0.3fMs/s).\n',sample_rate/1e6);
end
% Downsample received signal
nSamples = ceil(ofdmInfo.SamplingRate/round(sample_rate)*size(Data,1));
nRxAnts = size(Data,2);
Data_downsample = zeros(nSamples, nRxAnts);
for i=1:nRxAnts
    Data_downsample(:,i) = resample(Data(:,1),ofdmInfo.SamplingRate, round(sample_rate));
end
%% Cell Search
fprintf('\nStart cell search...\n');

duplexModes = {'TDD' 'FDD'}; %Currently build for FDD signal but keep TDD cell search.
cyclicPrefixes = {'Normal' 'Extended'}; %Currently build for Normal CP but keep extended CP situation for future.

searchcell.MaxCellCount = 2;
searchcell.SSSDetection = 'PostFFT';
peakMax = -Inf;

% Here we might see mulitple Cell IDs in the list, we will display all of
% them but only keep the first one for next step;
for duplexMode = duplexModes
    for cyclicPrefix = cyclicPrefixes
        enb.DuplexMode = duplexMode{1};
        enb.CyclicPrefix = cyclicPrefix{1};
        [enb.NCellID, offset, peak] = lteCellSearch(enb, Data_downsample, searchcell);
        offset = offset(1);
        peak = peak(1);
        if(peak > peakMax)
            enbMax = enb;
            offsetMax = offset;
            peakMax = peak;
        end
    end
end

fprintf('This is %s signal, All the cell we saw from signal:', enbMax.DuplexMode);
disp(enbMax.NCellID(1));

enb = enbMax;
offset = offsetMax;
enb.NCellID = enb.NCellID(1);
%% Show the confidence of cell search result
corr = cell(1,3);
for i = 0:2
    enb.NCellID = floor(enb.NCellID/3)*3 + mod(enb.NCellID+i,3);
    [~,corr{i+1}] = lteDLFrameOffset(enb, Data_downsample); %Use PSS and SSS to check frame offset.
    corr{i+1} = sum(corr{i+1},2);
end
threshold = 1.3 * max([corr{2}, corr{3}]); %1.3 is empirically obtained.
if (max(corr{1}) < threshold)
    warning('Synchronization signal correlation was weak; Cell ID might be wrong');
end
enb.NCellID = enbMax.NCellID(1);
enb.NSubframe = 0;

%% CFO correction
fprintf('\nCorrect Frequency Offset...\n');
Data_downsample = Data_downsample(offset+1:end,:);
if (strcmpi(enb.DuplexMode,'TDD')) %Not dealing with TDD yet
    enb.TDDConfig = 0;
    enb.SSC = 0;
end
delta_f = lteFrequencyOffset(enb, Data_downsample);
fprintf('Frequency offset: %0.3fHz\n', delta_f);
Data_downsample = lteFrequencyCorrect(enb, Data_downsample, delta_f);

%% OFDM demodulation 
fprintf('\nOFDM demodulation...\n');
enb.CellRefP = 4; % Assume 4 cell-specific reference signals initially for decoding attempt. Will change in the future.
resourceGrid = lteDLResourceGridSize(enb); %NRB * Nsymbol * Antport
Nsymbol = resourceGrid(2); %Number of symbolsl in one subframe
rxgrid = lteOFDMDemodulate(enb, Data_downsample);
if(isempty(rxgrid))
    fprintf('Signal not correct, no subframe successfully be decoded\n');
    return;
end

%% Channel Estimation
cec.PilotAverage = 'UserDefined';
cec.FreqWindow = 13;
cec.TimeWindow = 9;
cec.InterpType = 'cubic';
cec.InterpWindow = 'Centered';
cec.InterpWinSize = 1;
[hest, nest] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:Nsymbol,:)); %hest is channel estimation and nest is noise estimation.

%% PBCH decoding, MIB parsing.
% MIB is carried on the BCH. This is transmitted with a fixed coding and
% modulation scheme and can be decoded after the initial cell search
% procedure. MIB directs UE to read SIB information which contains
% information for basic cell access parameters.
fprintf('\nDecode BCH...\n');
pbchIndices = ltePBCHIndices(enb);
[pbchRxSym, pbchHestSym] = lteExtractResources(pbchIndices, rxgrid(:,1:Nsymbol,:), hest(:,1:Nsymbol,:,:));
[pbchBits, pbchSym, pbchNFmod4, mib, cellrefP] = ltePBCHDecode(enb, pbchRxSym, pbchHestSym, nest);
enb.CellRefP = cellrefP;
enb = lteMIB(mib, enb);

enb.NFrame = pbchNFmod4;
disp(enb);

if(enb.CellRefP==0 || enb.NDLRB==0)
    fprintf('MIB decoding failed\n');
    return;
end
%% OFDM demodulate with new enb
fprintf('\nOFDM demodulation with new enb...\n');
ofdmInfo = lteOFDMInfo(enb);
if (sample_rate~=ofdmInfo.SamplingRate)
    if (sample_rate < ofdmInfo.SamplingRate)
        warning('The received signal sampling rate (%0.3fMs/s) is lower than the desired sampling rate for cell search / RNTI decoding (%0.3fMs/s); cell search / RNTI decoding may fail.',sameple_rate/1e6,ofdmInfo.SamplingRate/1e6);
    end
    fprintf('\nResampling from %0.3fMs/s to %0.3fMs/s for cell search / RNTI decoding...\n',sample_rate/1e6,ofdmInfo.SamplingRate/1e6);
else
    fprintf('\nResampling not required; received signal is at desired sampling rate for cell search / RNTI decoding (%0.3fMs/s).\n',sample_rate/1e6);
end
nSamples = ceil(ofdmInfo.SamplingRate/round(sample_rate)*size(Data,1));
Data_resampled = zeros(nSamples, nRxAnts);
for i = 1:nRxAnts
    Data_resampled(:,i) = resample(Data(:,i), ofdmInfo.SamplingRate, round(sample_rate));
end

delta_f = lteFrequencyOffset(enb, Data_resampled);
fprintf('Frequency offset: %0.3fHz\n', delta_f);
Data_resampled = lteFrequencyCorrect(enb, Data_resampled, delta_f);

%Find beginning of frame
offset = lteDLFrameOffset(enb, Data_resampled);
Data_resampled = Data_resampled(1+offset:end,:);

resourceGrid = lteDLResourceGridSize(enb); %NRB * Nsymbol * Antport
Nsymbol = resourceGrid(2); %Number of symbolsl in one subframe
rxgrid = lteOFDMDemodulate(enb, Data_resampled);
if(isempty(rxgrid))
    fprintf('Signal not correct, no subframe successfully be decoded\n');
    return;
end
%% Decode RNTI information from SIB 
% The SIB1 is in subframe 5 of every even frame, Few steps to get to decode RNTI
% 1. Demodulate PCFICH which carries CFI. CFI will indicate how many
% subframes is used for PDCCH allocation.
% 2. PDCCH decoding.
% 3. Blind RNTI search in PDCCH.
while(size(rxgrid,2) > 0)
    pbchIndices = ltePBCHIndices(enb);
    [pbchRxSym, pbchHestSym] = lteExtractResources(pbchIndices, rxgrid(:,1:Nsymbol,:), hest(:,1:Nsymbol,:,:));
    [pbchBits, pbchSym, pbchNFmod4, mib, cellrefP] = ltePBCHDecode(enb, pbchRxSym, pbchHestSym, nest);
    enb.CellRefP = cellrefP;
    enb = lteMIB(mib, enb);
%     
    enb.NFrame = pbchNFmod4;
    
    if(enb.CellRefP==0 || enb.NDLRB==0)
        fprintf('MIB decoding failed\n');
        return;
    end
    enb.NSubframe = 0;
    rnti_dict_subframe = dictionary;
    rxframe = rxgrid(:,1:Nsymbol*10,:);
    rxframe_save = rxgrid(:,1:Nsymbol*10,:);
    hest_mat = [];
    nest_mat = [];
    while(size(rxframe,2) > 0)
        rxsubframe = rxframe(:,1+Nsymbol*enb.NSubframe:Nsymbol*(enb.NSubframe+1),:);
        [hest, nest] = lteDLChannelEstimate(enb, cec, rxsubframe);
        hest_mat = [hest_mat,hest];
        nest_mat = [nest_mat,nest];
        pcfichIndices = ltePCFICHIndices(enb);
        [pcfichRxSym, pcfichHestSym] = lteExtractResources(pcfichIndices, rxsubframe, hest);
        pcfichBits = ltePCFICHDecode(enb, pcfichRxSym, pcfichHestSym, nest);
        cfi = lteCFIDecode(pcfichBits);
        enb.CFI = cfi;
        pdcchIndices = ltePDCCHIndices(enb);
        [pdcchRxSym, pdcchHestSym] = lteExtractResources(pdcchIndices, rxsubframe, hest);
        [pdcchBits, pdcchSym] = ltePDCCHDecode(enb, pdcchRxSym, pdcchHestSym, nest);
        ueConfig.ControlChannelType='PDCCH';
        ueConfig.EnableCarrierIndication = 'Off';
        ueConfig.SearchSpace = 'UESpecific';
        ueConfig.EnableMultipleCSIRequest = 'off';
        ueConfig.EnableSRSRequest = 'off';
        [RNTI_dict, candidate_Rnti, Active_Rnti] = dci_rnti_decode(enb, ueConfig, pdcchBits, candidate_Rnti, Active_Rnti);
        rnti_info = struct;
        rnti_info.rawrnti = RNTI_dict;
        rnti_info.pdcchbits = pdcchBits;
        rnti_info.resolve = false;
        rnti_info.rnti = dictionary;
        rnti_dict_subframe(enb.NSubframe) = rnti_info;
        if(mod(enb.NSubframe,10) == 9)
            for i = enb.NSubframe - 9: enb.NSubframe
                [conflict_flag, potential_rnti, conflict_rnti] = ActiveRnti_filter(rnti_dict_subframe(i), Active_Rnti);
                if (conflict_flag == false && numEntries(potential_rnti) > 0)
                    rnti_dict_subframe(i).resolve = true;
                    rnti_dict_subframe(i).rnti = potential_rnti;
                    continue;
                end
                if (conflict_flag == false && numEntries(potential_rnti) == 0)
                    rnti_dict_subframe(i).resolve = false;
                    continue;
                end

                if (conflict_flag == true)
                    [conflict_flag, potential_rnti, conflict_rnti] = pdcch_filter(enb, rxframe_save(:,1+Nsymbol*i:Nsymbol*(i+1),:), potential_rnti, hest_mat(:,1+Nsymbol*i:Nsymbol*(i+1),:), nest_mat(i+1), conflict_flag, conflict_rnti);
                end
                if (conflict_flag == false && numEntries(potential_rnti) > 0)
                    rnti_dict_subframe(i).resolve = true;
                    rnti_dict_subframe(i).rnti = potential_rnti;
                    continue;
                end
                if (conflict_flag == false && numEntries(potential_rnti) == 0)
                    rnti_dict_subframe(i).resolve = false;
                    continue;
                end

                if (conflict_flag == true)
                    [conflict_flag, potential_rnti] = Active_freq_filter(conflict_rnti, Active_Rnti, potential_rnti);
                end
                if (conflict_flag == false && numEntries(potential_rnti) > 0)
                    rnti_dict_subframe(i).resolve = true;
                    rnti_dict_subframe(i).rnti = potential_rnti;
                    continue;
                end
                if (conflict_flag == false && numEntries(potential_rnti) == 0)
                    rnti_dict_subframe(i).resolve = false;
                    continue;
                end

%                 if (~isKey(rnti_subframe,i))
%                     tmp_rnti = dictionary;
%                     rnti_sf_dict = rnti_dict_subframe(i);
%                     if (rnti_sf_dict.numEntries > 0)
%                         rnti_list = keys(rnti_sf_dict);
%                         for j = 1:length(rnti_list)
%                             if(iscached(appeared_rnti,num2str(rnti_list(j))))
%                                 tmp_rnti(rnti_list(j)) = rnti_sf_dict(rnti_list(j)).Allocation;
%                             end
%                         end
%                         rnti_subframe(i) = tmp_rnti;
%                     end
%                 end
            end
            break;
        end
%         rxframe(:,1:Nsymbol,:) = [];
        enb.NSubframe = mod(enb.NSubframe + 1,10);
    end
    for i = enb.NSubframe - 9: enb.NSubframe
        if rnti_dict_subframe(i).resolve == false
            fprintf('%d\t0\t0\n\n',enb.NFrame*10+i);
            continue;
        end
        rnti_subframe = rnti_dict_subframe(i).rnti;
        rnti_list = keys(rnti_subframe);
        for m = 1:length(rnti_list)
            fprintf('%d\t%d\t',enb.NFrame*10+i,rnti_list(m));
            disp(rnti_subframe(rnti_list(m)).dcimsg.Allocation);
        end
    end
    
    enb.NSubframe = mod(enb.NSubframe + 1,10);
    rxgrid(:,1:Nsymbol*10,:) = [];
end






