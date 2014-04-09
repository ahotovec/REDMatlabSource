function [varargout] = pickDates2arraypicks(year,day,varargin)
% Takes PICKDATES .mat files according to input year, day, channel, etc...
% and constructs the set of array-wide observations, according to user 
% input critiera concerning differential observation times, number of 
% recievers observing an event, and time window. No normalization is done
% on the seismograms. No filtering is performed.
%
% USAGES
% [S,T] = pickDates2arraypicks(year,day,opt)
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
% S:    The coral structure of seismograms corresponding to array-wide
%       observations
% T:    The array wide observation date vectors
%
% EXAMPLE: Find the coral structure and date vectors corresponding to
% seismic events observed on all 6 stations of an array on Julian day 346,
% day 2006, using .EPZ pick date files. Extract 10 sec of data for each
% observation. Allow 1 second differential travel times between events:
%
% opt.cutTime = 5; 
% opt.level = 6-1; 
% opt.rsec = 1;
% opt.chan = EPZ;
% opt.cutTime = 5;
% opt.level = 6;
% opt.rsec = 1;
%
% [S,T]=pickDates2arraypicks(2005,346,opt);
%
%-----------------------------------------------------------------------
% Latest Edit: 25.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

%set default values for opt structure
opt.level       = 5;
opt.cutTime     = 5;
opt.rsec        = 1;
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

%get the cell of p-wave pick times in datecell. Note that getPickDates
%ensures day is a 3-digit number, say a user input of 1 makes 001, so that
%getPickDates reads in the correct file.
[datecell]  = getPickDates(day, year, opt.chan);

%Find where datecell has garbage entries and discard them.
%This is faster than the vectorized form for a 6 station array.  Better to
%use than cellfun.
for k=1:length(datecell)
    datecut                         = datecell{k};
    datecut(:,datenum(datecut')==0) = []; %datenum only used here to identify zeros
    datecell{k}                     = datecut;
end;

%get the picks that show up on level number of stations at the same time
C = clusterdates(datecell,opt.rsec,opt.level);

%Determine the number of individual observations per level
numobs =  cellfun('size', C, 2);

%discard empty values in staCell and C.
C=C(numobs>0);

%Retain only non-empty observations exclusive to a fixed number of stations
C = uniquepickDates(C);
I = cellfun(@isempty,C);

%Get pick times for the array-wide observed events in order to output for
%use with a picker.
if(~isempty(C)),
    C=C(~I); 
    T=C{end};
end;

if(isempty(C)),
    T = [];
    S = [];
end;

%Make seismogram matrix for array-wide observed events
%cutTime: amount of seismogram to keep
%---------------------------------------------------------------------

Sref        = setCoralFields;
allFields   = fieldnames(Sref);
allFields   = setdiff(allFields,'staCode');

for k=1:size(T,2)
    if(k==1)
        [S,Sacs]    = pascSac2Coral('SAC','SAC',timeadd(T(:,1),-opt.cutTime),timeadd(T(:,1),+opt.cutTime));
        S           = coralSetNames(S,Sacs);
        S           = coralTaper(coralDetrend(S));
    else
        [temp,Sacs] = pascSac2Coral('SAC','SAC',timeadd(T(:,k),-opt.cutTime),timeadd(T(:,k),+opt.cutTime) );
        temp        = coralSetNames(temp,Sacs);
        temp        = coralTaper(coralDetrend(temp));
        S           = cat(2,S,temp);
    end;
end;

S = coralSetRefs(S,Sref,'staCode',allFields{:});

%% Save arraypicks to a .mat file for clustering, etc... later
arraypicks.pTime        = T;
arraypicks.coralstruc   = S;
saveFile                = sprintf('%s.%s.%i.%.3i.mat','ARRAYPICKS',opt.chan,year,day);
disp(saveFile) 
save(saveFile,'arraypicks');

if(nargout>=1), varargout{1} = T;           end;
if(nargout>=2), varargout{2} = S;           end;
if(nargout>=3), varargout{3} = saveFile;    end;
return;