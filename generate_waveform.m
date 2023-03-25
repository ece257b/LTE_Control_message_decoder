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
cfg.PDSCH.RNTI = 255;
cfg.PDSCH.Rho = 0;
cfg.PDSCH.RVSeq = [0 1 2 3];
cfg.PDSCH.NHARQProcesses = 8;
cfg.PDSCH.PMISet = 1;
cfg = lteRMCDL(cfg);

% input bit source:
in = [1; 0; 0; 1];


% Generation
[waveform, grid, cfg] = newlteRMCDLTool(cfg, in);
repeatwaveform = waveform;
Fs = cfg.SamplingRate; 								 % Specify the sample rate of the waveform in Hz

fiq = fopen('~/Downloads/lte_test_4000.dat','w');
fwrite(fiq, waveform, 'float');
fclose(fiq);

for i=1:4000
    cfg = struct('RC', 'R.2', ...
        'DuplexMode', 'FDD', ...
        'NCellID', 4, ...
        'TotSubframes', 10, ...
        'NumCodewords', 1, ...
        'Windowing', 0, ...
        'AntennaPort', 1, ...
        'NSubframe', i*10);
    
    cfg.OCNGPDSCHEnable = 'Off';
    cfg.OCNGPDCCHEnable = 'Off';
    cfg.PDSCH.TxScheme = 'Port0';
    cfg.PDSCH.RNTI = 255;
    cfg.PDSCH.Rho = 0;
    cfg.PDSCH.RVSeq = [0 1 2 3];
    cfg.PDSCH.NHARQProcesses = 8;
    cfg.PDSCH.PMISet = 1;
    cfg = lteRMCDL(cfg);
    
    % input bit source:
    in = [1; 0; 0; 1];
    
    if mod(i,100) == 0
        fprintf("%d\n",i);
    end
%     
%     % Generation
    [waveform, grid, cfg] = newlteRMCDLTool(cfg, in);
    Fs = cfg.SamplingRate; 								 % Specify the sample rate of the waveform in Hz

    fiq = fopen('~/Downloads/lte_test_4000.dat','a');
    fwrite(fiq, waveform, 'float');
    fclose(fiq);
end