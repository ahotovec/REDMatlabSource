%Build an organized catalog out of the results of the template search
clear all;close all;clc
%size of MasterStack
sizel=[1 444];

cd /Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog

%combine each year first and save in catalog folder
if 0
for year=2003:2013
    Dall=cell(sizel);
    evtimes=cell(sizel);
    A=dir(['../RCMtemplsearch/',num2str(year),'*']);
    tvecT=[]; total=[]; %initialize things
    switchy=1; %switch=1 if want to combine multiplets, otherwise 0 if just want to get multperhr
    for i=1:length(A)
        load(['../RCMtemplsearch/',A(i).name]);
        [a, b]=multPerHr(D,1);
        total=[total, b];
        tvecT=[tvecT, a];
        if switchy==1
            for j=1:length(D)
                temp=D{j};
                if ~isempty(temp)
                    xc=[temp.xcorval];
                    temp=temp(isfinite(xc));
                    Dall{j}=[Dall{j} temp];
                    evtimes{j}=[evtimes{j}; datenum([temp(:).recStartTime]')];
                end
            end
        end
    end
    [tvecT, I]=sort(tvecT);
    total=total(I);
    save(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/Yearly/',num2str(year),'Dall.mat'],'Dall','tvecT','total','evtimes','-v7.3')
end
end
%%
%initialize empty event mat files
if 0
Drcm=[]; evtimes1=[]; numev=[];starttime=[];endtime=[];name={};stack=[];
for i=1:max(sizel)
    save(['IndividualSets/',num2str(i),'.mat'],'Drcm','evtimes1','numev','starttime','endtime','name','stack')
end
end
%% Then go through and stuff each set into matlab structures separately
if 0
year=2003:2013;
tvecAll=[];
totalAll=[];
for i=1:length(year)
    load(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/Yearly/',num2str(year(i)),'Dall.mat'])
    for j=1:length(Dall)
        if ~isempty(Dall{j})
            load(['IndividualSets/',num2str(j),'.mat'])
            %add in data with data already there for this event
            Drcm=[Drcm; Dall{j}'];
            evtimes1=[evtimes1; evtimes{j}];
            numev=length(Drcm);
            starttime=min(evtimes1);
            endtime=max(evtimes1);
            name=datestr(starttime,'yymm-ddHHMMSS');
            %resave event
            save(['IndividualSets/',num2str(j),'.mat'],'Drcm','evtimes1','numev','starttime','endtime','name','stack')
        end
    end
    %combine tvecT and total with all others
    tvecAll=[tvecAll tvecT];
    totalAll=[totalAll total];
end
[tvecAll, I]=sort(tvecAll);
totalAll=totalAll(I);
save('multperhrAll.mat','tvecAll','totalAll');
end

%% Then go through each event and make stacks and compute summary vectors for things
if 1
%initiatilize catalog summary vectors
numevAll=zeros(sizel);
numev7=zeros(sizel);
names=cell(sizel);
domfreqstack7=zeros(sizel);
startdates=zeros(sizel);
enddates=zeros(sizel);
signalwidth7=zeros(sizel);
kurtosF7=zeros(sizel);
kurtosT7=zeros(sizel);
meanRecint=zeros(sizel);
stdRecint=zeros(sizel);
meanRecint7=zeros(sizel);
stdRecint7=zeros(sizel);
freqVar=zeros(sizel); %variance of frequency content (second moment)

A=dir('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/*.mat');

for i=1:length(A)
    load(['IndividualSets/',A(i).name]);
    num=str2num(A(i).name(1:end-4)); %get number of event
    %extract event specific info
    xcorvals=[Drcm.xcorval]; %xcorvals relative to stack
    recint=[0;diff(sort(evtimes1))].*86400; %recurrence interval in seconds
    recint7=[0; diff(sort(evtimes1(xcorvals>=0.7)))].*86400;
    domfreq=coralDominantFreq(Drcm,0);
    numev=length(evtimes1);
    starttime=datestr(min(evtimes1));
    endtime=datestr(max(evtimes1));
    %stack of all and spectrum of each stack
    Drcm=coralPad(Drcm);
    stackAll=coralRaStack(Drcm');
    [f, a]=amplitudespectrum(coralTvec(stackAll),stackAll.data);
    
    %stack of those with 0.7 and up (if there aren't any, delete this set
    %cause it sucks
    if max(xcorvals>=0.7)
        stack7=coralRaStack(Drcm(xcorvals>=0.7)');
        [f7, a7]=amplitudespectrum(coralTvec(stack7),stack7.data);
    else
        delete(['IndividualSets/',A(i).name])
        continue
    end
    
    %stack of those with 0.8 and up
    if max(xcorvals>=0.8)
        stack8=coralRaStack(Drcm(xcorvals>=0.8)');
        [f8, a8]=amplitudespectrum(coralTvec(stack8),stack8.data);
    else
        stack8=[];
        f8=[];a8=[];
    end
    
    %save individual info
    save(['IndividualSets/',num2str(num),'.mat'],'Drcm','evtimes1','numev','starttime','endtime','stackAll','recint','recint7','domfreq','f','a','stack7','f7','a7','stack8','f8','a8')
      
    %extract event summary info
    numevAll(num)=length(evtimes1);
    numev7(num)=length(evtimes1(xcorvals>=0.7));
    names{num}=datestr(min(evtimes1),'yymm-ddHHMMSS');
    startdates(num)=min(evtimes1);
    enddates(num)=max(evtimes1);
    meanRecint(num)=mean(recint);
    stdRecint(num)=std(recint);
    meanRecint7(num)=mean(recint7);
    stdRecint7(num)=std(recint7);
    evtimes{num}=evtimes1;
    domfreqstack7(num)=coralDominantFreq(stack7,0);
    signalwidth7(num)=coralSignalWidth(stack7);
    kurtosF7(num)=kurtosis(a7)-3;
    kurtosT7(num)=kurtosis(stack7.data)-3; %normalized kurtosis

    %freqVar(num)=ones(sizel); %variance of frequency content (second moment)

end

%save summary info
save('CatalogSummary.mat','numevAll','numev7','names',...
    'startdates','enddates','meanRecint','stdRecint','meanRecint7',...
    'stdRecint7','domfreqstack7','signalwidth7','kurtosF7','kurtosT7','evtimes');
end
%% Make set of stacks
%Stack7s=zeros(sizel);
if 1
evtimes=cell(sizel);

A=dir('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/*.mat');
for i=1:length(A)
    load(['IndividualSets/',A(i).name],'stack7','evtimes1');
    num=str2num(A(i).name(1:end-4)); %get number of event
    Stack7s(num)=stack7;
    evtimes{num}=evtimes1;
end
save('Stack7s.mat','Stack7s')
save('evtimes.mat','evtimes')
end