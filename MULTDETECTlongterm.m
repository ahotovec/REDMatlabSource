% Script to run multiplet analysis using Josh Carmichael's RaCorrelation codes
clear all
close all
clc
%add path to all the scripts required (edit to match their location on your
%own computer

addpath /Users/Kate/MATLAB1/Rainier_all/LongTermSearch
addpath /Users/Kate/MATLAB1/Kates/MultCodes/MultipletPrograms
%javaaddpath /Users/Kate/MATLAB1/Kates/usgs.jar
javaaddpath /Users/Kate/MATLAB1/Kates/IRIS-WS-latest.jar

home='/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/';
cd(home);

%loads parameter structure to be used in following MATLAB functions
%Edit and double check defaultOptStruc before running blindly
opt = defaultOptStruc;

%enter name of the stations to use in the analysis in the form
%stations={'STAR.EHZ.UW'}; 
stations={'RCM.EHZ.UW'};
sta=stations{1}(1:end-7);


%enter starting and ending month and year
startYear=2012;
startMonth=11; %number
endYear=2013;
endMonth=6; %stops at the 1st of this month
 
%make vector of starts of each month in datenum format
monthvec=monthvector(startYear,startMonth,endYear,endMonth);
%make a directory to save results for this station, will just print a
%warning that it already exists if it already exists
mkdir([home,sta,'multsearch']);

% if 1
%matlabpool(3)
for i=1:length(monthvec)-1
    mon=monthvec(i);
%change back to home directory if not already there
cd([home,sta,'multsearch'])
    
%make a new directory for the month and year you are running
mkdir(datestr(mon,'yymm'))
%change to the new directory
cd(datestr(mon,'yymm'))
    
%enter starttime and endtime for the dates you want to search over, use
%at least half days
%format: startTime=[YYYY MM DD 00 00 00];
startTime=datevec(monthvec(i));
endTime=datevec(monthvec(i+1));

% Gets dates in the right formats
year=startTime(1);
H=datenum([startTime(1) 01 01 00 00 00]); %get day number of first day of year
startdayjulian=round(datenum(startTime)-H+1); %get julian day of start
enddayjulian=round(datenum(endTime)-H+1); %get julian day of end
dayvec=startdayjulian:enddayjulian; %make a vector of the julian days you will run the program over for file naming purposes

% This next program loops over each day. It gathers the seismograms from the winston server for the
%events corresponding to all of your pick times from the previous step. It
%determines which of these "picks" show up on enough stations within a short
%enough time interval (you define these in defaultOptStruc) to be called an
%event. It saves only the events for all stations in
%ARRAYPICKS.COMP.YYYY.DDD.mat,
%with one file for each day.

getPicks2(stations,startTime,endTime,opt);

end
%matlabpool close
%end
%%
%matlabpool(2)
for i=1:length(monthvec)-1
    mon=monthvec(i);
%change back to home directory if not already there
cd([home,sta,'multsearch'])
    
%change to the new directory
cd(datestr(mon,'yymm'))
    
%enter starttime and endtime for the dates you want to search over, use
%at least half days
%format: startTime=[YYYY MM DD 00 00 00];
startTime=datevec(monthvec(i));
endTime=datevec(monthvec(i+1));

% Gets dates in the right formats
year=startTime(1);
H=datenum([startTime(1) 01 01 00 00 00]); %get day number of first day of year
startdayjulian=round(datenum(startTime)-H+1); %get julian day of start
enddayjulian=round(datenum(endTime)-H+1); %get julian day of end
dayvec1=startdayjulian:enddayjulian; %make a vector of the julian days you will run the program over for file naming purposes

%make clusters for each day
for m = dayvec1(1:end-1)
    try
        junk=dir(['CLUSTERPICKS*',num2str(year),'.',num2str(m,'%.3i'),'.mat']);
        if isempty(junk)
[S,T] = getDailyClusters(year,m,stations,opt); 
        else
          fprintf([junk(1).name,' already exists\n'])  
        end
    catch
        fprintf(['Failed to make CLUSTERPICKS file for ',num2str(m,'%.3i'),' year ',num2str(startTime(1)),'. that stinks\n'])
        continue
    end
end;

%
%This step takes all the multiplets detected each day and compares them
%with each other, then it takes stacks of big sets, makes that a template,
%and pulls out any straggler multiplet events. It saves results as
%MULTTEMPL*.mat, one file for the entire experiement period.
%multtempl.multiplets contains the data for all multiplets grouped into
%clusters, multtempl.template contains the stack of each set to use as a
%template to find more
%try
    junk=dir(['MULTTEMPL*',num2str(year),'*.mat']);
            if isempty(junk)
    [multtempl]=getRaClusters(opt);
            else  
              fprintf([junk(1).name,' already exists\n'])  
            end
% catch
%     fprintf(['Failed to make MULTTEMPL file for month ',num2str(startTime(2)),' year ',num2str(startTime(1)),'. that sucks\n'])
%     
%     continue
% end

 end
 %matlabpool close



