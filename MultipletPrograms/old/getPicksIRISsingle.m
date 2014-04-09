function getPicksIRISsingle(stationnames,startTime,endTime,opt)
% function to find picks and save the data for a single station
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
    
    %for m=1:length(stationnames)
        
        recStart= datevec(round(H+day-1)); 
        date1   = recStart';
        date2   = timeadd(date1,timeWin);
        
        pTime   = {};
        
        for n = 1:(opt.nWins);
            
             S     = getIRISdata(stationnames,date1',date2');      
             S     = coralDemean(S);             
            if opt.lpfilt>0
                S      = coralTaper(S);
                S=coralFilter(S,opt.lpfilt,'low',2,'zero');
            end

            [picks] = pickerRatio(coralDetrend(S), opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, opt.test);          
            
            pTime   = cat(1,pTime,picks);     
            
                for j=1:size(picks,1)
                    optcut.absStartTime=(timeadd(picks{j},-opt.cutTimeB4));
                    optcut.absEndTime=(timeadd(picks{j},opt.cutTimeAF));
                    
                    temp = coralCut(S,optcut);
                    temp = coralTaper(coralDetrend(temp));
                    if j==1 && n==1; 
                        D=temp;
                    else
                        D=cat(2,D,temp);
                    end
                end;
                
            if( n==(opt.nWins) )
                               
                arraypicks.pTime        = pTime;
                arraypicks.coralstruc   = D;
                eventcount=length(pTime);
                if day>365; day=day-365; end
                saveFile                = sprintf('%s.%s.%i.%.3i.mat','ARRAYPICKS',opt.chan,date1(1),day);
                disp(saveFile)
                save(saveFile,'arraypicks','eventcount');

            else
                
            date1   = date2;
            date2   = timeadd(date1,timeWin);
            end
            
        end
        
end

disp('Done');
return;


