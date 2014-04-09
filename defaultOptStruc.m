function opt = defaultOptStruc
% This script simply outputs a default opt structure for usage in functions
% like raPicks2raCluster, pickDates2arraypicks, etc... There are no input
% arguments, and like setCoralFields, it exists in the raw data directory,
% as opposed to in SeisTools, as each raw data directory will have a
% separate default opt structure.  It can be edited for convenience, and is
% most useful to expediate using an opt structure with a large set of
% diverse fields.
opt.datatype = 'iris'; %'iris' or 'winston'

opt.winston     = 'mazama'; %name of the winston server to download data from if using winston

opt.thresh      = 0.7; %correlation threshold value for clustering, use lower value for multiple stations
opt.maxNumClust = 200;  %maximum number of clusters computed without using
                        %loops (memory constraint)
                        
%opt.IndXCorr=0.7; %x-corr threshold for individual station x-correlations
opt.localTime   = 8*(3600); %hours to add to UTC time to get local time at array %changes from 7 to 8 on Nov7,2010 but I don't care
opt.cutTimeB4     = 5; %the number of seconds to cut before a date/time pick
opt.cutTimeAF    = 20; %number of seconds to cut after a date/time pick

opt.lpfilt      = 0;%15;%10; %Hz %Low pass filter cutoff, must be zero for no lowpass
opt.bpfilt     =[1 10];%,must be zero for no bandpass
opt.chan        = 'EHZ'; %the channel used to compute picks

opt.level       = 1;    %the number of stations required to see an event for clustering/array picking, doesn't matter if you only have one station

opt.rsec        = 3; %the maximum number of seconds between p-wave arrivals within 
                       %the array before declaring an event.
opt.tbin        = 3600; %the number of seconds events counts are binned
opt.nWins       = 4; %daily windows used to pick events (prevents memory issues)
opt.sWindow     = 0.8; %short-time window used for p-wave picking
opt.lWindow     = 7; %long time winow used for p-wave picking
opt.rCutoff     = 3; %STA/LTA cut off
opt.nInterval   = 10; %the minimum time interval (in secs) required between distinct events
opt.taumax    =  15;   %maximum signal width allowed, see coralSignalWidth for details (recommended 15 sec for typical multiplets)
opt.test        = 0; %boolean value set to false
