function getPicks2(stationnames,startTime,endTime,opt)
% function to find picks and save the data for multiple stations and
% calculate hourly rms for each station. This version is different from
% getPicks in that it cuts out events with a signal width greater than that
% specified in opt.taumax (uses coralSignalWidth.m to do this)
%
% USAGE
% getPicks2(stationnames,startTime,endTime,opt)
%modified from Josh Carmichael's (josh.carmichael@gmail.com)
%codes by Kate Allstadt for use with data accessed from
%IRIS DMC or from Winston waveserver and combining pickdates 
%and arraypicks steps
%
% INPUT
% stationnames: station name in form {'RCS.EHZ.UW'}
% startTime: start vector [2010 01 01 00 00 00]
% endTime: end vector
% opt:  structure of input options specified in defaultOptStruc.
%
% opt.chan: channel used
% opt.sWindow:      short-time window over which to pick using STA/LTA
% opt.lWindow:      long-time window over which to pick using STA/LTA
% opt.rCutoff:      short time to long term window ratio cut off;
%               
% nInterval:    Number of seconds between events.  Wes recommends 10.
% opt.nWins:        The number of windows to cut each Sac file into prior to
%               use with pickerRatio.  Default is 6 for 6, 4 hours windows
% opt.taumax

%OUTPUTS
%FILL IN HERE
%-----------------------------------------------------------------------

H=datenum([startTime(1) 01 01 00 00 00]);
startdayjulian=round(datenum(startTime)-H+1);
enddayjulian=round(datenum(endTime)-H+1);
numSec  = 24*60*60; %number of seconds in a day
timeWin = numSec/(opt.nWins); %number of seconds in each time window
hrrms.datevec=[];
hrrms.rms=[];

dayvec=startdayjulian:enddayjulian;
%dayvec(dayvec>365)=dayvec(dayvec>365)-365;
dayvec=dayvec(1:end-1);
optcut.cutType='absTime';

for k=1:length(dayvec)
    %try
    listing = dir(['ARRAYPICKS*',num2str(startTime(1)),'*',num2str(dayvec(k),'%.3i'),'*']);
    if ~isempty(listing)  %if ARRAYPICKS file already exists for this day, skip to next day
        fprintf(['Arraypicks file for julian day',num2str(dayvec(k)),' year ',num2str(startTime(1)),' already exists, skipping to next day\n'])
        continue
    else
    
    day     = dayvec(k);
    sprintf('Picking day %i',char(day))
    
    recStart= datevec(round(H+day-1));
    date1   = recStart';
    date2   = timeadd(date1,timeWin);
    
    pTime   = {};
    
    for n = 1:(opt.nWins);
     
        %get data for all stations
        if strcmpi(opt.datatype,'winston')
            S     = coralWinData(stationnames,date1',date2',opt.winston);
        elseif strcmpi(opt.datatype,'iris')
            S     = getIRISdata(stationnames,date1',date2');
        end
        S     = coralDemean(S);
        Sorig=S; %save original data
        [rmstemp datevectemp]=hourlyrms(S);
        
        if length(opt.bpfilt)==2
            S      = coralTaper(S);
            S=coralFilter(S,opt.bpfilt,'bandpass',2,'minimum');
        elseif opt.lpfilt>0 && opt.bpfilt==0
            S=coralTaper(S);
            S=coralFilter(S,opt.lpfilt,'low',4,'minimum');
        end
        
        if length(stationnames)>1
        for m=1:length(stationnames)
            [picks] = pickerRatio(coralDetrend(S(m)), opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, opt.test);
            
            datecell{1,1,m}   = [picks{:}];
        end
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
                
                temp = coralCut(Sorig,optcut)'; %save original unfiltered data
                temp = coralTaper(coralDetrend(temp));
                if j==1 && n==1;
                    D=temp;
                else
                    D=cat(2,D,temp);
                end
            end;
            
            if n==1
                pTime=datenum(T');
            else
                pTime=[pTime; datenum(T')];
            end
            
        elseif length(stationnames)==1
            [picks] = pickerRatio(coralDetrend(S), opt.sWindow, opt.lWindow, opt.rCutoff, opt.nInterval, opt.test);
            pTime   = cat(1,pTime,picks);
            
            if length(picks)~=1 %only do this step if there are pick times
            for j=1:size(picks,1)
                optcut.absStartTime=(timeadd(picks{j},-opt.cutTimeB4));
                optcut.absEndTime=(timeadd(picks{j},opt.cutTimeAF));
                temp = coralCut(Sorig,optcut); %save original data
                temp = coralTaper(coralDetrend(temp));
                if j==1 && n==1;
                    D=temp;
                else
                    D=cat(2,D,temp);
                end
            end;
            end
        end
        if isempty(pTime)
            D=[];
        end
        
        %calculate hourly rms for current data
        hrrms.rms=[hrrms.rms; rmstemp];
        hrrms.datevec=[hrrms.datevec datevectemp];
        clear rmstemp datevectemp
        
        %get rid of signals with too long duration
        %only do this step if D exists
        if exist('D')
        if ~isempty(D)
        tau=coralSignalWidth(D);
        I=find(tau<15);
        D=D(I);
        pTime=pTime(I);
        end
        else 
            continue
        end
        
        if( n==(opt.nWins) ) 
            eventcount=length(pTime);
            
            if eventcount<500;
            arraypicks.pTime        = pTime;
            arraypicks.coralstruc   = D;
            
            if day>365; dayp=day-365; else dayp=day;end
            saveFile                = sprintf('%s.%s.%i.%.3i.1.mat','ARRAYPICKS',opt.chan,date1(1),dayp);
            disp(saveFile)
            save(saveFile,'arraypicks','eventcount','hrrms');   
                        
            else %divide arraypicks into two files so clusterpicks doesn't get overwhelmed

                %first file
                halfpt=round(length(pTime)/2);
                arraypicks.pTime=pTime(1:halfpt);
                arraypicks.coralstruc=D(:,1:halfpt);
                eventcount=length(arraypicks.pTime);
                if day>365; dayp=day-365; else dayp=day;end
                saveFile                = sprintf('%s.%s.%i.%.3i.1.mat','ARRAYPICKS',opt.chan,date1(1),dayp);
                disp(saveFile)
                save(saveFile,'arraypicks','eventcount','hrrms');
                %second file
                arraypicks.pTime=pTime(halfpt:end);
                arraypicks.coralstruc=D(:,halfpt:end);
                eventcount=length(arraypicks.pTime);
                saveFile                = sprintf('%s.%s.%i.%.3i.2.mat','ARRAYPICKS',opt.chan,date1(1),dayp);
                disp(saveFile)
                save(saveFile,'arraypicks','eventcount','hrrms'); 
                
            end
            
        else   
            date1   = date2;
            date2   = timeadd(date1,timeWin);
        end      
    end
    end
%     catch
%         fprintf(['failed to get arraypicks file for julian day ',num2str(dayvec(k)),' year ',num2str(startTime(1)),'\n'])
%         continue
%     end
end
disp('Done');
return;


