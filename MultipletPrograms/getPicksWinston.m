function getPicksWinston(stationnames,startTime,endTime,opt)
% function to find picks and save the data for multiple stations
% USAGE
% getPicksIRISsingle(stationnames,startTime,endTime,opt)

%modified from Josh Carmichael's (josh.carmichael@gmail.com) 
%clustering codes by Kate Allstadt for use with data accessed from 
%IRIS DMC
%
% INPUT
% stationnames: station name in form {'RCS.EHZ.UW'}
% startTime: start vector [2010 01 01 00 00 00]
% endTime: end vector
% opt.chan:  the channel for the SAC file
%
% 
% opt.sWindow:      short-time window over which to pick using STA/LTA
% opt.lWindow:      long-time window over which to pick using STA/LTA
% opt.rCutoff:      short time to long term window ratio cut off; 
%               Wes reccomends 2.
% nInterval:    Number of seconds between events.  Wes recommends 10.
% opt.nWins:        The number of windows to cut each Sac file into prior to
%               use with pickerRatio.  Default is 6 for 6, 4 hours windows
%-----------------------------------------------------------------------

H=datenum([startTime(1) 01 01 00 00 00]);
startdayjulian=round(datenum(startTime)-H+1);
enddayjulian=round(datenum(endTime)-H+1);
numSec  = 24*60*60; %number of seconds in a day
timeWin = numSec/(opt.nWins); %number of seconds in each time window
         
dayvec=startdayjulian:enddayjulian;
%dayvec(dayvec>365)=dayvec(dayvec>365)-365;
dayvec=dayvec(1:end-1);
optcut.cutType='absTime';

for k=1:length(dayvec)    
    day     = dayvec(k);
    sprintf('Picking day %i',char(day))

        recStart= datevec(round(H+day-1)); 
        date1   = recStart';
        date2   = timeadd(date1,timeWin);
        
        pTime   = {};
        
        for n = 1:(opt.nWins);
            
                
                %get data for all stations
                S     = coralWinData(stationnames,date1',date2',opt.winston);
                S     = coralDemean(S);
                if length(opt.bpfilt)==2
                    S      = coralTaper(S);
                    S=coralFilter(S,opt.bpfilt,'bandpass',2,'minimum');
                elseif opt.lptfilt>0 && opt.bpfilt==0
                    S=coralTaper(S);
                    S=coralFilter(S,opt.lpfilt,'low',4,'minimum');
                end
                
                for m=1:length(stationnames)
                [picks] = pickerRatio(coralDetrend(S(m)), opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, opt.test);
                
                
                datecell{1,1,m}   = [picks{:}];
                end
                %pTime   = cat(1,pTime,picks);
                
                %once all stations are picked through, find the ones that are "arraypicks"
                
                %get the picks that show up on level number of stations at the same time
                C = clusterdatesK(datecell,opt.rsec,opt.level);
                clear datecell
                
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
                    T=C{:};
                end;
                
                if(isempty(C)),
                    T = [];
                    S = [];
                end;
                    
                for j=1:size(T,2)
                    optcut.absStartTime=(timeadd(T(:,j),-opt.cutTimeB4));
                    optcut.absEndTime=(timeadd(T(:,j),opt.cutTimeAF));
                    
                    temp = coralCut(S,optcut)';
                    temp = coralTaper(coralDetrend(temp));
                    if j==1 && n==1;
                        D=temp;
                    else
                        D=cat(2,D,temp);
                    end
                end;
                
                if n==1; pTime=datenum(T');
                else
                pTime=[pTime; datenum(T')];
                end
                %end
                
                if( n==(opt.nWins) )
                    
                    arraypicks.pTime        = pTime;
                    arraypicks.coralstruc   = D;
                    if day>365; day=day-365; end
                    saveFile                = sprintf('%s.%s.%i.%.3i.mat','ARRAYPICKS',opt.chan,date1(1),day);
                    disp(saveFile)
                    save(saveFile,'arraypicks');
                    
                else
                    
                    date1   = date2;
                    date2   = timeadd(date1,timeWin);
                end
                
            
            
        end
end
disp('Done');
return;


