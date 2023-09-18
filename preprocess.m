theFilename = ['Subj002_Sess001_run001.snirf']; %Change this line to process a different file
acquired = SnirfClass(['F:\paper\real_data\snirf_converted_raw\' 'Subj002_Sess001_run001.snirf']);
tmpacquired = SnirfClass();
%The real "data" matrix is in acquired.data.dataTimeSeries
%And display some information about it.
tmpacquired.stim = copy(acquired.stim);

for i = 1:size(tmpacquired.stim, 2)   
    % Clear other attributes in the stimulus group
    tmpacquired.stim(1,i).data = []; % Clear the data attribute
    tmpacquired.stim(1,i).name = []; % Clear the name attribute.
    tmpacquired.stim(1,i).states = [];
end

startIndex = 4;  % Starting Index
endIndex = 45;    % End Index

% Create a new array of structs containing only the elements to be retained
newStim = [];
for temp0 = 1:numel(tmpacquired.stim)
    if temp0 < startIndex || temp0 > endIndex
        newStim = [newStim, tmpacquired.stim(temp0)];
    end
end

% Updating the stim property in tmpacquired
tmpacquired.stim = newStim;

for temp1 = 1:5
  onset = acquired.stim(1,2*temp1 - 1).data(1);
  duration = acquired.stim(1,2*temp1).data(1) - acquired.stim(1,2*temp1 - 1).data(1);
  tmpacquired.stim(1,1).name = '1';
  tmpacquired.stim(1,1).data(temp1,1:3) = [onset,duration, 1];
  tmpacquired.stim(1,1).states(temp1,1:2) = [onset,duration];
end

for temp2 = 1:5
  onset = acquired.stim(1,8 + 3*temp2).data(1);
  duration = acquired.stim(1,10 + 3*temp2).data(1) - acquired.stim(1,8 + 3*temp2).data(1);
  tmpacquired.stim(1,2).name = '2';
  tmpacquired.stim(1,2).data(temp2,1:3) = [onset,duration, 1];
  tmpacquired.stim(1,2).states(temp2,1:2) = [onset,duration];
end

for temp3 = 1:5
  onset = acquired.stim(1,22 + 4*temp3).data(1);
  duration = acquired.stim(1,25 + 4*temp3).data(1) - acquired.stim(1,22 + 4*temp3).data(1);
  tmpacquired.stim(1,3).name = '3';
  tmpacquired.stim(1,3).data(temp3,1:3) = [onset,duration, 1];
  tmpacquired.stim(1,3).states(temp3,1:2) = [onset,duration];
end
acquired.stim = copy(tmpacquired.stim);
acquired.Info();


% Get the size of measurementlist
[row, col] = size(acquired.data.measurementList);

% Loop over the first row of the measurementlist
for i = 1:col
    acquired.data.measurementList(1, i).sourcePower = 0;  % Assign sourcePower to 0
    acquired.data.measurementList(1, i).moduleIndex = i;  % Assign moduleIndex to n.
end
numSources = size(acquired.probe.sourcePos3D, 1);  % Get Number of sources

acquired.probe.sourceLabels = cell(numSources, 1);  % Create an empty array of cells matching the number of sources
for iSrc = 1:numSources
    acquired.probe.sourceLabels{iSrc} = sprintf('S%d', iSrc);
end
numDetectors = size(acquired.probe.detectorPos3D, 1);  % Get the number of detectors

acquired.probe.detectorLabels = cell(numDetectors, 1);  % Create an empty array of cells matching the number of detectors

for iDet = 1:numDetectors
    acquired.probe.detectorLabels{iDet} = sprintf('D%d', iDet);
end

plotOptions.shortChannelDistance = 1.5; %as used below in the GLM for regressing out
plotOptions.stim = acquired.stim;%Short channel, parameter 1.5

baselineDuration = acquired.stim(1, 1).data(1); % Duration of baseline period (time before the first stimulus)
baselineEndTime = baselineDuration; % End time of baseline period

numChannels = size(acquired.data.dataTimeSeries, 2); % Get the number of channels

% Initialize variables for min and max signal
minSignal = inf;
maxSignal = -inf;
% Initialize variable for average SNR
minSNR = inf;
maxSNR = -inf;

for iCh = 1:numChannels
    % Get the time series data for the current channel
    timeSeries = acquired.data.dataTimeSeries(:, iCh, :);
    
    % Find the index corresponding to the end time of the baseline period
    baselineEndIdx = find(acquired.data.time <= baselineEndTime, 1, 'last');
    
    % Extract the baseline data for the current channel
    baselineData_channel = timeSeries(1:baselineEndIdx, :);
    
    % Calculate the average for each wavelength
    baselineMean = mean(baselineData_channel, 1);
    
    % Update min and max signal
    minSignal = min(minSignal, min(baselineMean));
    maxSignal = max(maxSignal, max(baselineMean));
    
    % Calculate SNR for each wavelength
    baselineStd = std(baselineData_channel, 1);
    SNR_wavelengths = baselineMean ./ baselineStd;
    
    % Update min and max SNR
    minSNR = min(minSNR, min(SNR_wavelengths));
    maxSNR = max(maxSNR,max(SNR_wavelengths));
end
param.SNRthresh = minSNR;
% Set param.dRange based on minSignal and maxSignal
param.dRange = [minSignal maxSignal];
param.mlActMan = [];
param.tIncMan = [];
param.SDrange = [0, 45.0];
param.mlActAuto = hmrR_PruneChannels(acquired.data,acquired.probe, ...
    param.mlActMan,param.tIncMan,param.dRange, param.SNRthresh, param.SDrange);

plotOptions.blockActive = param.mlActAuto;
%Assigning values to logical arrays
 [hfig,hBGAxis,hChAxis] = myHomer3_plotSnirfDataOnProbe(acquired.data, acquired.probe,plotOptions);
 
 outResults.dod = hmrR_Intensity2OD(acquired.data);  %Convert light intensities into optical densities
disp(' == Plotting: 04_OD_RawOpticalDensities')

%[hfig,hBGAxis,hChAxis] = myHomer3_plotSnirfDataOnProbe(outResults.dod, acquired.probe,plotOptions);

outResults.dod = hmrR_BandpassFilt(outResults.dod, 0.005, 0.35);

%[hfig,hBGAxis,hChAxis] = myHomer3_plotSnirfDataOnProbe(outResults.dod, acquired.probe,plotOptions);

param.turnon = 1;
outResults.dod = hmrR_PreprocessOD_LinearFit(outResults.dod,param.turnon);
%disp(' == Plotting: 06_OD_Detrended')


param.tMotion = 0.5;
param.tMask = 2;
param.STDEVthresh = 14;
param.AMPthresh =  maxSignal - minSignal ;
[param.tIncAuto,param.tIncAutoCh] = hmrR_MotionArtifactByChannel(outResults.dod,acquired.probe, ...
    param.mlActMan,param.mlActAuto,param.tIncMan,param.tMotion,param.tMask,param.STDEVthresh,param.AMPthresh);
disp(' == Plotting: 07_OD_MotionArtefactDetection1')

%param.trange = [-20 20];%此处参数需要修改
%[ExcludedStim,~] = hmrR_StimRejection(outResults.dod,acquired.stim,...
    %param.tIncAuto,param.tIncMan,param.trange);

%disp(' == Plotting: 08_OD_StimRejection')

%plotOptions.stim = ExcludedStim;
%plotOptions.blockActive = param.mlActAuto;
tIncAuto_1 = sum(param.tIncAuto{1,1} == 1);
tIncAuto_0 = sum(param.tIncAuto{1,1} == 0);
ratio = tIncAuto_0 / (tIncAuto_1 + tIncAuto_0);

param.p=0.99;
param.turnon = 1; %Watch out! Overwriting the above turnon
outResults.dod = hmrR_MotionCorrectSpline (outResults.dod, ...
    param.mlActAuto,param.tIncAutoCh,param.p,param.turnon);

disp(' == Plotting: 09_OD_MotionCorrection')

[param.tIncAuto,param.tIncAutoCh] = hmrR_MotionArtifactByChannel(outResults.dod,acquired.probe, ...
    param.mlActMan,param.mlActAuto,param.tIncMan,param.tMotion,param.tMask,param.STDEVthresh,param.AMPthresh);
    % This method can work on optical densities or on
    % concentrations. Here I'm applying it over optical
    % densities.

disp(' == Plotting: 10_OD_MotionArtefactDetection2')
tIncAuto_1 = sum(param.tIncAuto{1,1} == 1);
tIncAuto_0 = sum(param.tIncAuto{1,1} == 0);
ratio1 = tIncAuto_0 / (tIncAuto_1 + tIncAuto_0);

outResults.dc  = hmrR_OD2Conc(outResults.dod, acquired.probe, [6.25,6.25]);
disp(' == Plotting: 11_Conc_RawMBLLReconstruction')



%outResults.dcAvg = hmrR_BlockAvg (outResults.dc, acquired.stim, [acquired.stim(1, 1).states(1, 1),acquired.stim(1, 1).states(5, 1) + acquired.stim(1, 1).states(5, 2)]);

[hfig,hBGAxis,hChAxis] = myHomer3_plotSnirfDataOnProbe(outResults.dc, acquired.probe,plotOptions);

%outResults.dcAvg = hmrR_BlockAvg (outResults.dc, acquired.stim, param.trange);
%disp(' == Plotting: 12_Conc_BlockAveraged')

%[hfig,hBGAxis,hChAxis] = myHomer3_plotSnirfDataOnProbe(outResults.dcAvg, acquired.probe,plotOptions);