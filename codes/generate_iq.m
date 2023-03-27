RNTI = 255;
subframeNumb=20;
filename='../lte_test_R2_2frame.dat';
cfg = struct('RC', 'R.2', ...
    'DuplexMode', 'FDD', ...
    'NCellID', 4, ...
    'TotSubframes', 10, ...
    'NumCodewords', 1, ...
    'Windowing', 0, ...
    'AntennaPort', 1, ...
    'NSubframe', 0);

cfg.OCNGPDSCHEnable = 'Off';
cfg.OCNGPDCCHEnable = 'Off';
cfg.PDSCH.TxScheme = 'Port0';
cfg.PDSCH.RNTI = RNTI;
cfg.PDSCH.Rho = 0;
cfg.PDSCH.RVSeq = [0 1 2 3];
cfg.PDSCH.NHARQProcesses = 8;
cfg.PDSCH.PMISet = 1;
cfg = lteRMCDL(cfg);

% input bit source:
in = [1; 0; 0; 1];


% Generation
[waveform, grid, cfg] = newlteRMCDLTool(cfg, in);
Fs = cfg.SamplingRate; 								 % Specify the sample rate of the waveform in Hz

wiq = zeros(size(waveform,1)*2,1);
wiq(1:2:end) = real(waveform(:,1));
wiq(2:2:end) = imag(waveform(:,1));
fiq = fopen(filename,'w');
fwrite(fiq, wiq, 'float');
fclose(fiq);

for i=1:subframeNumb-1
    cfg = struct('RC', 'R.2', ...
        'DuplexMode', 'FDD', ...
        'NCellID', 4, ...
        'TotSubframes', 10, ...
        'NumCodewords', 1, ...
        'Windowing', 0, ...
        'AntennaPort', 1, ...
        'NSubframe', 0, ...
        'NFrame', i);
    
    cfg.OCNGPDSCHEnable = 'Off';
    cfg.OCNGPDCCHEnable = 'Off';
    cfg.PDSCH.TxScheme = 'Port0';
    cfg.PDSCH.RNTI = RNTI;
    cfg.PDSCH.Rho = 0;
    cfg.PDSCH.RVSeq = [0 1 2 3];
    cfg.PDSCH.NHARQProcesses = 8;
    cfg.PDSCH.PMISet = 1;
    cfg = lteRMCDL(cfg);
    
    % input bit source:
    in = [1; 0; 0; 1];
    
%     if mod(i,100) == 0
%         fprintf("%d\n",i);
%     end
%     
%     % Generation
    [waveform, grid, cfg] = newlteRMCDLTool(cfg, in);
    Fs = cfg.SamplingRate; 								 % Specify the sample rate of the waveform in Hz
    wiq = zeros(size(waveform,1)*2,1);
    wiq(1:2:end) = real(waveform(:,1));
    wiq(2:2:end) = imag(waveform(:,1));
    fiq = fopen(filename,'a');
    fwrite(fiq, wiq, 'float');
    fclose(fiq);
end