function opt = redOptions

% set all the options in opt structure for RED

% General options

opt.timeFile = 'timeFile.txt'; %name of text file that says when the last run was started (should be in path used for outputs in redMAIN unless you specify the full file path here)
opt.archiveTime = 60; %time, in number of days, until a family expires and is sent into archive (waveforms deleted)..this is measured in time since the last event


%% Options for event picking
opt.timelimit = 12/24;  %time, in days, of maximum amount of data that can be loaded in for picking at a time
opt.cutB4 = 3; %number of seconds to cut before pick time
opt.cutAF = 12; %number of seconds to cut after pick time
opt.sWindow = 0.8; %      short-time window over which to pick using STA/LTA
opt.lWindow = 10; %      long-time window over which to pick using STA/LTA
opt.rCutoff = 3; %     short term to long term window ratio cut off; 
opt.dataType = 'IRIS'; %IRIS or winston
opt.winston = 'products01'; %winston waveserver to use if opt.dataType == 'winston'
opt.nInterval = opt.cutAF; %Number of seconds required between events.  Wes recommends 10

opt.Kwindow = [1 15]; %time window over which to compute kurtosis
opt.Kmin = 3.5; %minimum kurtosis allowed
opt.Kmax = 200; %maximum kurtosis allowed
opt.ORlim = 0.06; %threshold for ratio of outliers (to delete bad events), 0.05 recommended

% filtering options for picking events
opt.useFilter = 1; %0 for no filter, 1 for yes filter
opt.filterType = 'bandpass'; %use coralFilter conventions
opt.cutoffFreq = [1 10]; %use coral convention
opt.order = 2;
opt.phase = 'minimum'; 

%% Options for event clustering

opt.dsamp = 20; %new samplerate (in Hz) for downsampling (to save on computing time), lp filtering performed to avoid aliasing

opt.mlag = 3; %+-lag time in seconds for cross correlation



