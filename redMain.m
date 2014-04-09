function redMain(station,optfile,outloc,timeflag,timeperiod)
%redMain(station,optfile,outloc,timeflag,timeperiod)
%INPUTS
%   station: String of station name to pick. Just one station. Must be in form {'STA.CHAN.NET'} (e.g. {'RCS.EHZ.UW'})
%   optfile: path (string) to location and name of file containing all the
%               options you defined for this search
%   outloc: path (string) to the directory where output files should be
%               saved (should probably be specific to your station)
%   timeflag:   1 if you will specify the start and end times
%               2 if you will specify the end time and period before the end
%                   time to scan
%               3 if you want to use now as the end time and specify the time
%                   before present to scan
%               4 if you want to use now as the end time and search for a file
%                   that specifies the last time the program was run that will be
%                   used to define the start time
%   timeperiod: the time period over which you want to search for repeaters
%               if flag=1, [starttime endtime] in any matlab date format
%               if flag=2, [endtime days] days specifies how far
%                   back you want to go in time in # of days
%               if flag=3, just input the number of days back you
%                   want to scan
%               if flag=4, just put a 0, doesn't matter. The name
%                  of the file used to save the time is specified in the redOptions.m file
%               

%OUTPUTS
%   no direct outputs from this function, but many files are output by the codes run from the
%   main script
%EXPLAIN THE FILES
%
%RED codes written by Kate Allstadt (allstadt.k@gmail.com)
%clustering algorithms based on codes written by Josh Carmichael

%% Initialization
%read in the options file
opt = optfile;

%change to the directory where the files will be located
cd(outloc)

% %put station name into a cell
% station={station};

%initialize folders needed (if they don't already exist)
if ~isdir('orphanTank') %this is where the orphan events will be kept
    mkdir('orphanTank')
end
if ~isdir('familyTank') %this is where event families will be stored
    mkdir('familyTank')
end
% if ~isdir('temp') %this is where the event pick files will be saved temporarily, then deleted once events are farmed out to orphan and family tanks
%     mkdir('temp')
% end
if ~isdir('figures') %this is where the magic happens
    mkdir('figures')
end
if ~isdir('archive') %this is where old families go to die once they have not had any new events for the specified amount of time (waveforms not saved?? Just metadata)
    mkdir('archive')
end

%get the start and end time
if timeflag==1  %start and end time were specified
    
    starttime = datenum(timeperiod(1));
    endtime = datenum(timeperiod(2));
    
elseif timeflag==2   %end time and period before end time were specified
    
    endtime = datenum(timeperiod(1));
    starttime = endtime-timeperiod(2);
    
    
elseif timeflag==3  %use now as end time and search back specified number of days
    endtime=now;
    starttime=now-timeperiod;
elseif timeflag==4    %use now as end time and go back till the last time the analysis was run (last end time used)
    %search for opt.timeFile
    A=dir(opt.timeFile);
    if isempty(A) || length(A)>1
        error('%s not found, or more than one file found',opt.timeFile)
    end
    %read in time file
    fid=fopen(opt.timeFile);
    starttime=fscanf(fid,'%f');
    fclose(fid);
    endtime=now; %define end time
    %new timeFile with this ending time will be written at the end of this
    %function if the run was successful
else
    error('You must define a timeflag between 1 and 4, fool. Try again')    
end


%% 
%find event pick times in the specified period and extract data
D = redPicker(station,datevec(starttime),datevec(endtime),opt);

% %save picks in temp folder
% save([outloc,'/temp/picks_',datestr(starttime,'yyyymmdd_HH:MM'),'-',datestr(endtime,'yyyymmdd_HH:MM'),'.mat'],'D')

%compare these events to existing families in familyTank (first to stacks
%with lower xcorr coeff, then if passes that test, xcorr against all other
%ones to see if it fits, and recluster into subevents etc.



%compare remaining events to each other and create new families, move into
%familyTank
%filter and downsample
samprate=1./D(1).recSampInt;
newrate=samprate./opt.dsamp;
Dlp=coralFilter(D,ceil(newrate./2),'low'); %LP filter to half of new sample rate
options.sintr_target = 1./newrate;
Dds=coralResample(Dlp,options);

[~, ~, ~, ~, maxC] = coralCrossCorr(Dds, opt); %(skip old step of coralCluster, manually input maxC to xcorMutClust)
[Scluster,Sorphan]    = xcorMutClust(S,opt,maxC);

%compare new events with orphans in orphanTank and pull any newly created
%families into familyTank



%move all unmatched events into orphanTank



%move any non active families into the archives



%delete any poor expired orphans from orphanTank



%update opt.timeFile with new endtime since scan was successful
fid=fopen(opt.timeFile,'w');
fprintf(fid,'%f',endtime);
fclose(fid);


end