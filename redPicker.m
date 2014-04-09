function D = redPicker(station,starttime,endtime,opt)

% function to find picks and save the data for multiple stations and
% calculate hourly rms for each station
%
% USAGE
% D = redPicker(station,starttime,endtime,opt)
%
% modified from Josh Carmichael's (josh.carmichael@gmail.com)
% clustering codes by Kate Allstadt for use with data accessed from
% IRIS DMC or from Winston waveserver and combining pickdates
% and arraypicks steps and for compatibility with RED
%
% INPUT
% station:          station name (just one) in form {'RCS.EHZ.UW'}
% starttime:        start time in any matlab date form
% endTime:          end time in any matlab date form
% opt:              structure of input options specified in redOptions.
%
% opt.sWindow:      short-time window over which to pick using STA/LTA
% opt.lWindow:      long-time window over which to pick using STA/LTA
% opt.rCutoff:      short time to long term window ratio cut off;
%                   Wes reccomends 2.
% opt.nInterval:    Number of seconds between events.  Wes recommends 10.
% opt.dataType:     'winston' to get data from winston (need to set opt.winston also if this is chosen), or 'IRIS' to get data
%                   from IRIS
% opt.timelimit:    Time maximum duration of data that can be loaded in
%                   (this is to save on data transfer time, can get bogged down if too long)
% opt.cutB4:        Number of seconds to cut out from data before picktime
%                   (also size of buffer)
% opt.cutAF:        Number of seconds to cut out from data after picktime
%                   (also size of buffer)
% opt.Kwindow:      window length for computing kurtosis
% opt.Kmin:         minimum kurtosis allowed
% opt.Kmax:         maximum kurtosis allowed
% opt.ORlim:        maximum percentage of outliers allowed (to eliminate
%                   square waveforms)
%
% OUTPUT
% D:                structure array containing data cut out around picktimes

%% Initialize
numSec  = 24*60*60; %number of seconds in a day
%hrrms.tvec=[]; %NOT SURE IF WILL KEEP THIS
%hrrms.rms=[];
optcut.cutType = 'absTime';
endtime=datenum(endtime);
starttime=datenum(starttime);
duration = endtime-starttime; %duration of analysis in days

%determine number of windows for downloading data
if opt.timelimit<duration %determine if more than one download window will be required
    nWin=ceil(duration/opt.timelimit);
else
    nWin=1; %only one time window needed
end


%% let's get some picktimes!

for n = 1:nWin;
    
    %determine start and end times and add a buffer of opt.cutB4 and
    %opt.cutAF to ends so data doesn't get cut off if a pick is right at
    %the end (no picks will be allowed to originate in this buffer though)
    if n==nWin %either end or there is only one nWin
        if n==1; %this is when there is only one window
            date1=datevec(starttime-opt.cutB4./numSec);
        else %this is when it's the last window
            date1=datevec(datenum(date2)-opt.cutB4./numSec-opt.cutAF./numSec);
        end
        date2=datevec(endtime+opt.cutAF./numSec);
    elseif n==1 && nWin>1 %it's the first window
        date1=datevec(starttime-opt.cutB4./numSec);
        date2=datevec(starttime+opt.timelimit+opt.cutAF./numSec);%don't need to add buffer because it takes old buffer already from last date2
    else %window in between
        date1=datevec(datenum(date2)-opt.cutB4./numSec-opt.cutAF./numSec); %take date2 from last window as start time of next (after subtracting buffer)
        date2=datevec(datenum(date2)+opt.timelimit);%don't need to add buffer because it takes old buffer already from last date2
        
    end
    
    %get data
    if strcmpi(opt.dataType,'winston')
        S     = coralWinData(station,date1',date2',opt.winston);
    elseif strcmpi(opt.dataType,'iris')
        S     = getIRISdata(station,date1,date2);
    end
    %pre-process
    S = coralDemean(S);
    S = coralTaper(S);
    Sorig = S; %save original data to cut raw events out later
    %[rmstemp datevectemp]=hourlyrms(S);
    
    if opt.useFilter==1 %filter data for detecting if needed
        S = coralFilter(S,opt.cutoffFreq,opt.filterType,opt.order,opt.phase);
    end
    
    %     %calculate hourly rms for current data
    %     hrrms.rms=[hrrms.rms; rmstemp];
    %     hrrms.datevec=[hrrms.datevec datevectemp];
    %     clear rmstemp datevectemp
    
    %do picking
    [picks] = pickerRatio(S, opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, 0);
    %picks=datenum([picks{:}]');
    
    %cut out data
    for j=1:size(picks,1)
        %don't use picks that are in the buffer window
        if datenum(picks{j}')<datenum(date1)+opt.cutB4/numSec || datenum(picks{j}')>datenum(date2)-opt.cutAF/numSec
            continue
        end
        optcut.absStartTime=(timeadd(picks{j},-opt.cutB4));
        optcut.absEndTime=(timeadd(picks{j},opt.cutAF));
        temp = coralCut(Sorig,optcut); %save original data
        if j==1 && n==1;
            Dt=temp;
        else
            Dt=cat(2,Dt,temp);
        end
        
    end
    
    
end

%now eliminate junky picks based on kurtosis (get rid of spikes) and percentage of outliers in
%signal (get rid of calibration pulses and others with flat tops)

K=coralKurtosis(Dt,opt.Kwindow);
OR=coralOutlierRatio(Dt,0);

D=Dt(K>opt.Kmin & K<opt.Kmax & OR<opt.ORlim);
junk=Dt(K<=opt.Kmin | K>=opt.Kmax | OR>=opt.ORlim);


end