function [varargout]=sac2RaPicksWinston(stationnames,startTime,endTime,opt)

%edited 10/18/2011 by Kate Allstadt for use with data from a winston server


% Takes SAC data in current directory, and input year, channel, day vector
% and optional input structure to get pick dates from SAC data.
%
% USAGE
% [pTime] = sac2RaPicks(year,dayvec);
% [pTime] = sac2RaPicks(year,dayvec,opt);
% [pTime] = sac2RaPicks(year,dayvec,'STATION');
% [pTime] = sac2RaPicks(year,dayvec,nWins);
% [pTime] = sac2RaPicks(year,dayvec,opt,'STATION',nWins);
%
% INPUT
% opt.chan:       The channel for the SAC file
% year:     The year for the SAC file
% dayvec:   A vector of Julian days, of the form [120,122,155]
% STATION:  A character string giving the station name to pick over
%
% The below are fields of pickerRatio structure P. P is optional.  There
% are defaults that are input in the event no P is provided.
%
% sWindow:      short-time window over which to pick using STA/LTA
% lWindow:      long-time window over which to pick using STA/LTA
% rCutoff:      short time to long term window ratio cut off; 
%               Wes reccomends 2.
% nInterval:    Number of seconds between events.  Wes recommends 10.
% nWins:        The number of windows to cut each Sac file into prior to
%               use with pickerRatio.  Default is 6 for 6, 4 hours windows
% test:         A Boolean for plotting.
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 26.Jan.2010
% Replaced nWins with opt.nWins in case 3 inputs were supplied
%
% 12.April.2010
% Replaced ch with opt.chan since opt includes chan field (or should).

% 05.May.2010
% Replaced redundant getSacFiles with getFiles function.
%-----------------------------------------------------------------------

%Default values to input into pickerRatio.m if no opt structure is provided
% opt.nWins       = 6;
% opt.sWindow     = 0.5;
% opt.lWindow     = 2.5;
% opt.rCutoff     = 3.2;
% opt.nInterval   = 5;
% opt.test        = 0;
% opt.chan        = 'EPZ';
% staname         = 'SAC';
% pTime           = {};
% [Sref]          = setCoralFields;
% allFields       = fieldnames(Sref);
 % allFields       = setdiff(allFields,'staCode');
% 
% if(nargin>3)
%     
%     for k=1:length(varargin)
%         
%         if(isstruct(varargin{k}))
%             
%             opt         = varargin{k};
%             
%             if(isfield(opt,'nWins')),
%                 nWins       = opt.nWins;
%             end;
%             
%         end;
%         
%         if(isscalar(varargin{k}))
%             
%             nWins = varargin{k};
%             
%         end;
%         
%         if(ischar(varargin{k}))
%             
%             staname = varargin{k};
%             
%         end;
%         
%     end;
% end;

H=datenum([startTime(1) 01 01 00 00 00]);
startdayjulian=round(datenum(startTime)-H+1);
enddayjulian=round(datenum(endTime)-H+1);
year=startTime(1);

dayvec=startdayjulian:enddayjulian;
dayvec=dayvec(1:end-1);

for k=1:length(dayvec)
    
    day     = dayvec(k);
    sprintf('Picking day %i',char(day))
%    stationList = getFiles(staname,opt.chan,num2str(day),'SAC',num2str(year));

    for m=1:length(stationnames)
        
        temp        = stationnames(m);
        saveName    = textscan(temp{1},'%s','Delimiter','.');
        saveName    = saveName{1}(1);
        
        sprintf('Picking through station %s',char(stationnames(m)))
        
%        [S,nfiles]	= pascSac2Coral(char(saveName),opt.chan,-1,-1,day,year);
        
 %       sprintf('%i files found for %s, chan %s', length(nfiles), char(staname), opt.chan)
        
%         if( not(isstruct( S) ) ),
%             pTime = {};
%             saveFile    = sprintf('%s.%s.%s.%i.%.3i.mat','PICKDATES',char(saveName),char(opt.chan),year,day)
%             save(saveFile,'pTime');
%             continue;
%         end;
        
        numSec  = 24*60*60; %number of seconds in a day
        
        timeWin = numSec/(opt.nWins);
        
        if( timeWin <= opt.lWindow ),
            pTime = {};
            saveFile    = sprintf('%s.%s.%s.%i.%.3i.mat','PICKDATES',char(saveName),char(opt.chan),year,day)
            save(saveFile,'pTime');
            continue;
        end;
        
        recStart= datevec(round(H+day-1)); 
        date1   = recStart';
        date2   = timeadd(date1,timeWin);
        
        pTime   = {};
        
        for n = 1:(opt.nWins);
            
%             [S,Sacs]    = pascSac2Coral( char(saveName), opt.chan, date1, date2 );
%             [S]         = coralSetNames(S,Sacs);
%             [S]         = coralSetRefs(S,Sref,'staCode',allFields{:});

            [S]     = coralWindata(stationnames(m),date1',date2',opt.winston);

            [picks] = pickerRatio(coralDetrend(S), opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, opt.test);          
            
            pTime   = cat(1,pTime,picks);            
            
            if( n==(opt.nWins) )
                
                saveFile    = sprintf('%s.%s.%s.%i.%.3i.mat','PICKDATES',char(saveName),char(opt.chan),year,day)
                save(saveFile,'pTime');
                disp(sprintf('Save pTime, length %i, %s.%s, start %s, end %s', length(pTime), char(saveName), char(opt.chan), datestr(recStart',0), datestr( date2',0 )))
                pTime   = {};
                
            end;
            
            date1   = date2;
            date2   = timeadd(date1,timeWin);
            
        end;
        
    end;
    
end;

if(nargout==1)
    
    varargout{1}    = pTime;
    
end;
disp('Done');
return;