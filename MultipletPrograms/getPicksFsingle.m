function getPicksFsingle(stations,startTime,endTime,opt,Bproj)
% function to find picks and save the data for multiple stations and
% calculate hourly rms for a single station using F=detector. 
%
% USAGE
% getPicksFsingle(stationnames,startTime,endTime,opt)
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
dayvec=startdayjulian:enddayjulian;
dayvec=dayvec(1:end-1);
year=startTime(1);

for k=1:length(dayvec)
    %try
    listing = dir(['ARRAYPICKS*',num2str(startTime(1)),'*',num2str(dayvec(k),'%.3i'),'*']);
    if ~isempty(listing)  %if ARRAYPICKS file already exists for this day, skip to next day
        fprintf(['Arraypicks file for julian day',num2str(dayvec(k)),' year ',num2str(startTime(1)),' already exists, skipping to next day\n'])
        continue
    else
    
    day     = dayvec(k);
    sprintf('Picking day %i',char(day))
    
    [D]=Fpicker(datenum(startTime)+k-1,datenum(startTime)+k,stations,opt,Bproj);
    D=D{:};
    pTime=datenum(timeadd([D.recStartTime],-opt.cutTimeB4)');

    eventcount=length(pTime);
            
            if eventcount<500;
            arraypicks.pTime        = pTime;
            arraypicks.coralstruc   = D;
            
            if day>365; dayp=day-365; else dayp=day;end
            saveFile                = sprintf('%s.%s.%i.%.3i.1.mat','ARRAYPICKS',opt.chan,year,dayp);
            disp(saveFile)
            save(saveFile,'arraypicks','eventcount');   
                        
            else %divide arraypicks into two files so clusterpicks doesn't get overwhelmed

                %first file
                halfpt=round(length(pTime)/2);
                arraypicks.pTime=pTime(1:halfpt);
                arraypicks.coralstruc=D(:,1:halfpt);
                eventcount=length(arraypicks.pTime);
                if day>365; dayp=day-365; else dayp=day;end
                saveFile                = sprintf('%s.%s.%i.%.3i.1.mat','ARRAYPICKS',opt.chan,year,dayp);
                disp(saveFile)
                save(saveFile,'arraypicks','eventcount');
                %second file
                arraypicks.pTime=pTime(halfpt:end);
                arraypicks.coralstruc=D(:,halfpt:end);
                eventcount=length(arraypicks.pTime);
                saveFile                = sprintf('%s.%s.%i.%.3i.2.mat','ARRAYPICKS',opt.chan,year,dayp);
                disp(saveFile)
                save(saveFile,'arraypicks','eventcount'); 
                
            end
     
   
    end
%     catch
%         fprintf(['failed to get arraypicks file for julian day ',num2str(dayvec(k)),' year ',num2str(startTime(1)),'\n'])
%         continue
%     end
end
disp('Done');
return;


