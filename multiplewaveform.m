% Generated by MATLAB(R) 9.13 (R2022b) and LTE Toolbox 3.8 (R2022b).
% Generated on: 26-Feb-2023 12:13:58

%% Generating Downlink RMC waveform
% Downlink RMC configuration
FrameNum=2000;

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
len = length(repeatwaveform);
Fs = cfg.SamplingRate; 								 % Specify the sample rate of the waveform in Hz

cfg = struct('RC', 'R.2', ...
    'DuplexMode', 'FDD', ...
    'NCellID', 4, ...
    'TotSubframes', FrameNum*10, ...
    'NumCodewords', 1, ...
    'Windowing', 0, ...
    'AntennaPort', 1, ...
    'NSubframe', 10);
%     
cfg.OCNGPDSCHEnable = 'Off';
cfg.OCNGPDCCHEnable = 'Off';
cfg.PDSCH.TxScheme = 'Port0';
cfg.PDSCH.RNTI = 255;
cfg.PDSCH.Rho = 0;
cfg.PDSCH.RVSeq = [0 1 2 3];
cfg.PDSCH.NHARQProcesses = 8;
cfg.PDSCH.PMISet = 1;
cfg = lteRMCDL(cfg);
%     
%     % input bit source:
in = [0; 0; 0; 1; 1; 1];
% % Generation
[waveform, grid, cfg] = newlteRMCDLTool(cfg, in);
repeatwaveform=[repeatwaveform waveform];


%% Visualize
% Spectrum Analyzer
% spectrum = spectrumAnalyzer('SampleRate', Fs);
% spectrum(waveform);
% release(spectrum);


%% Transmit waveform over the air
masterClockRate = 15360000;
interpolationFactor = 1;
usrpTx = comm.SDRuTransmitter(Platform='B200');
usrpTx.SerialNum ='3272344';
usrpTx.CenterFrequency = 2000000000;
usrpTx.Gain = 89.75;
usrpTx.ChannelMapping = 1;
usrpTx.LocalOscillatorOffset = 1;
usrpTx.PPSSource = 'Internal';
usrpTx.ClockSource = 'Internal';
usrpTx.MasterClockRate = masterClockRate;
usrpTx.InterpolationFactor = interpolationFactor;
usrpTx.TransportDataType = 'int16';
usrpTx.EnableBurstMode = false;

% Transmit waveform (for 10 sec):
stopTime = 1000; % sec
t = 0; % sec
fprintf('Transmission started.\n');
while t<stopTime
    for i=1:FrameNum
        waveform = repmat(repeatwaveform(:,i),1,1);
        usrpTx(waveform);
        t = t+length(waveform)/Fs;
    end
end

fprintf('Transmission stopped.\n')
release(usrpTx);
