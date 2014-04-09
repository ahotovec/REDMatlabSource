function [varargout] = pickDates2countmetric(year,day,varargin)
% Takes PICKDATES .mat files according to input year, day, channel, etc...
% and constructs and bins a count score for array observations, according 
% to user input critiera concerning differential observation times, number 
% of recievers observing an event, and time window. Datenum is used
% internally, but with no consequence whatsoever on data interpretation.
%
% USAGES
% [binLoc]                  = pickDates2countmetric(year,day,opt);
% [binLoc,countBin]         = pickDates2countmetric(year,day,opt);
% [binLoc,countBin,fname]   = pickDates2countmetric(year,day,opt);
%
% INPUT
% year:     The four digit year corresponding to the pick dates.
% day:      The three digit Julian day corresponding to the pick dates.
%
% opt:      An optional input structure with the following fields:
% chan      The channel as it appears in the pick date save file/the coral
%           structure for the array.
% rsec      The maximum differential travel time between any station pair 
%           that indicates a single event was observed between that station 
%           pair
% level     The number stations used to declare an 'array-wide' 
%           observation
% cutTime   The half-width (in sec) of the time window used to store the 
%           array response function (coral structure)
%
% OUTPUT
% binLoc:   The temporal location of the bins
% dateBin:  The binned/histogrammed data
% fname:    The name of the file saved with fields
%
% %EXAMPLE: Find count score for all seismic events observed on all 6 
% %stations of an array on Julian day 346, day 2006, using .EPZ pick date 
% %files. Bin data by 1800 seconds (2 bins per hour). Allow 1 second 
% %differential travel times between events: cutTime = 5; level = 6; 
% %rsec = 1; The bin locations and count scores are output in T and S
% respectively.
% 
% opt.level       = 5;
% opt.cutTime     = 5;
% opt.rsec        = 1;
% opt.tbin        = 3600;
% opt.chan        = 'EPZ';
%
% [T,S]=pickDates2countmetric(2006,346,opt);
%
%-----------------------------------------------------------------------
% Latest Edit: 19.March.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 15.Nov.2009
% Authored function form of function
%
% 19.March.2010
% stationList2Cluster now has an opt.level option; datenum is still used,
% but only for binning; no stored data has been converted using datenum.
%
% 24.June.2010
% countmetric output now contains uniform output in field 'hist'. Empty
% values, if they exist, are retained to preserve structure of countmetric.
%-----------------------------------------------------------------------

%set default values for opt structure
opt.level       = 5;
opt.cutTime     = 5;
opt.rsec        = 1;
opt.tbin        = 3600;
opt.chan        = 'EPZ';

%intialize output
if(nargout > 1)
    for k = 1:nargout-1;
        varargout{k} = [];
    end;
end;

%check input
for k = 1:length(varargin)
    if(isstruct(varargin{k}))
        opt = varargin{k};
    end;
end;

%get the cell of p-wave pick times in datecell
[datecell,flist,stationList] = getPickDates(day,year,opt.chan);

%Find where datecell has garbage entries and discard them.
%This is faster than the vectorized form for a 6 station array.  Better to
%use than cellfun.
for k=1:length(datecell)
    datecut                         = datecell{k};
    datecut(:,datenum(datecut')==0) = []; %datenum is fine here--it's not replacing an actual date
    datecell{k}                     = datecut;
end;

%get the picks that show up on level number of stations at the same time
C           = clusterdates(datecell,opt.rsec,opt.level);

%make a corresponding organizational cell of station names
%staCell     = stationList2cluster(stationList);
staCell     = stationList2cluster(stationList,opt.level);

%Load blank coral structure containing geographical station coordinates
S           = setCoralFields;

%Get distance metric;
distmetric  = picks2distmetric(staCell,S);

%Determine the number of individual observations per level
%numobs      =  cellfun('size', C, 2);

%discard empty values in staCell and C.
%Nt=length(numobs(numobs>0));
%C       = C(numobs>0);
%staCell = staCell(numobs>0);
%distmetric = distmetric(numobs>0);

%Retain only non-empty observations exclusive to a fixed number of stations
C = uniquepickDates(C);
%I = cellfun(@isempty,C);

%Retain non-empty counts only:
%I           = not(logical(I));
%C           = C(I);
%distmetric    = distmetric(I); %now a matrix, not a cell
distmetric    = distmetric./max(distmetric);
%staCell     = staCell(I);

%get the date vector for the begining of the day
daystart    = [year;1;1;0;0;0];
daystart    = timeadd(daystart, (day-1)*24*3600 );

numBins     = floor((24*3600)./(opt.tbin));
tbin        = 24*3600/numBins;
binLoc      = timeadd(daystart, linspace(tbin/2,24*3600-tbin/2,numBins));
binLoc      = datenum( binLoc' );

%find the number of events for each time window
Cmetric     = cellfun(@transpose, C, 'UniformOutput', false);
Cmetric     = cellfun(@datenum, Cmetric, 'UniformOutput', false);

%% Compute the Count Metric
% This block of code bins the dates and multiplies the frequency by the
% distance metric.

countBin = zeros( size(binLoc,1), size(C,1) );
for k = 1:size(C,1)
    
    temp            = hist( Cmetric{k}, binLoc );
    temp            = temp.*[distmetric(k)];
    countBin(:,k)   = temp;
    
end;

%% Save countmetric to a .mat file
countmetric.binLoc      = binLoc;
countmetric.hist        = countBin;
countmetric.stations    = staCell;

saveFile = sprintf('%s.%s.%.3i.%s.%i.%.3i.mat','COUNTMETRIC',...
        'ARRAY',round((24*3600)/tbin),(opt.chan),year,day);

disp(saveFile) 
save(saveFile,'countmetric');

if(nargout>=1), varargout{1} = binLoc;      end;
if(nargout>=2), varargout{2} = countBin;    end;
if(nargout>=3), varargout{3} = saveFile;    end;

return;