function [MultAllSta]=getallsta(stalist,existmult,opt)
%This function retrieves data for the stations not used in detection listed
%in stalist. Gets raw unfiltered data for the same time intervals as those listed in
%existingmultiplets by loading in the data for the time period between
%earliesttime and latest time

MultAllSta=existmult;
numexiststa=size(existmult{1,1},1);

%find first and last multiplets and round down and up to get the time span
%to load in
for i=1:size(existmult,1)
    %get first and last start time of each multiplet
    dates1=[existmult{i,1}.recStartTime];
    dates1=datenum(dates1');
    mind(i)=min(dates1);
    maxd(i)=max(dates1);
end
startTime=floor(min(mind));
endTime=ceil(max(maxd));
%Load half a day at a time
NumHlfDays=(endTime-startTime)*2;

for i=1:NumHlfDays
    startt=datevec(startTime+(i-1)/2);
    endt=datevec(startTime+i/2);
    
    %get the data for all those stations
    if strcmpi(opt.datatype,'winston')
        D1     = coralWinData(stalist,startt,endt,opt.winston);
    elseif strcmpi(opt.datatype,'iris')
        D1     = getIRISdata(stalist,startt,endt);
    end

%     %calculate rms for all stations
%     [rmstemp datevectemp]=hourlyrms(D1);
    
    %Find empty data and fill with zeros
    %first, find a nonempty D1
    for n=1:size(D1,2)
        if ~isempty(D1(1,n).data)
            rst=D1(1,n).recStartTime;
            rnd=D1(1,n).recNumData;
            data=zeros(rnd,1);
            rsi=D1(1,n).recSampInt;
            break
        end
    end
    D1=coralPad(D1);
    for n=1:size(D1,2)
        if isempty(D1(1,n).recStartTime)
            D1(1,n).recStartTime=rst;
            D1(1,n).recNumData=rnd;
            D1(1,n).data=data;
            D1(1,n).recSampInt=rsi;
        end
    end
%     if opt.lpfilt>0
%         D1=coralDemean(D1);
%         D1=coralTaper(D1);
%         D1=coralFilter(D1,opt.lpfilt,'low',2,'zero');
%     end
    
    %Cut out the same pieces as in the existing multiplets
    for j=1:size(existmult,1)
        for k=1:size(existmult{j,1},2)
            
            opt1.absStartTime=existmult{j,1}(1,k).recStartTime;
            opt1.absEndTime=timeadd(existmult{j,1}(1,k).recStartTime,opt.cutTimeB4+opt.cutTimeAF);
            opt1.cutType='absTime';
            if datenum(opt1.absStartTime')>=datenum(startt) && datenum(opt1.absEndTime')<=datenum(endt)
                Dcut=coralCut(D1,opt1);
                Dcut=coralDemean(Dcut);
                %place in cell array
                for m=1:length(stalist)
                    %add blank xcorval field so can be combined with
                    %existing structures
                    Dcut(1,m).xcorval=[];
                    MultAllSta{j,1}(numexiststa+m,k)=Dcut(1,m);
                end
            end
        end
    end
%             %calculate hourly rms for current data
%         hrrms.rms=[hrrms.rms; rmstemp];
%         hrrms.datevec=[hrrms.datevec datevectemp];
%         clear rmstemp datevectemp
end

save('MultAllSta.mat','MultAllSta')%,'hrrms');
