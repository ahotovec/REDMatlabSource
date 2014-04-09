%reconcile the multdetect and template catalogs into one catalog including
%all of them

%written by Kate on 12 July 2013

%template catalog is king, add other detections missed by this method by
%adding any detections that aren't within timemin seconds of already existing
%detections and circshift those to line up with template ones

%load in full catalog from regular multiplet search (should be in same
%order as template events
if 0
clear all;clc;close all;
timemin=25; %max differnce in seconds require to declare separate events

load /Users/Kate/MATLAB1/Rainier_all/LongTermSearch/templsearchResults/MASTER_FULL_RCM_0.7thresh_1112-0513.mat
temp=MasterFull0p7thresh;
load /Users/Kate/MATLAB1/Rainier_all/LongTermSearch/templsearchResults/MASTER_FULL_RCM_0.7thresh.mat
%combine them
MasterFull=[MasterFull0p7thresh,temp];
clear MasterFull0p7thresh
MasterFull=MasterFull([1:789,791:end]);
sizel=length(MasterFull);
end

%%
%Go through each template event and add in missing events to template times
%(they should already be circshifted to the same times if they are the
%0.7thresh versions)
% if the template search didn't
%find any events for that set, fill in the events completely from the
%multiplet search and fill in the other info for the catalog in that case
if 0
for i=1:length(MasterFull)
    %find if there is a template set for this event set
    A=dir(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/',num2str(i),'.mat']);
    if isempty(A) %fill in event set completely from scratch
        Drcm=MasterFull{1,i}';
    else %just add detections to what is already there 
       %get event times from multdetect
       load(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/',A.name])
       if isempty(Drcm)
           continue
       end
       temp=MasterFull{1,i}';
       multDtimes=datenum([temp.recStartTime]');
       minT=zeros(size(multDtimes));
       for j=1:length(multDtimes) %find minimum time difference between each of these and all the template eventtimes
           minT(j)=min(abs(multDtimes(j)-evtimes1));
       end
       %keep any that have a time difference more than timemin 
       temp2=temp(minT>timemin./86400);
       Drcm=[Drcm; temp2];
       fprintf(['added ',num2str(length(temp2)),' out of ',num2str(length(MasterFull{1,i})),' events to ',num2str(i),'.mat\n']);
       %minT
    end
           save(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/',num2str(i),'.mat'],'Drcm')
end
end
%% update catalog 
% go through each event and make stacks and compute summary vectors for things
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

cd('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/')
A=dir('*.mat');

for i=1:length(A)
    load(A(i).name);
    if isempty(Drcm)
        continue
    end
    %line them all up
    opt.refseis=1;%only get cross corr relative to first event
    Drcm=coralPad(Drcm);
    [~, ~, ~, ccLag] = coralCrossCorr(Drcm, opt);
    %if median shift is greater than 5 seconds, it probably means the first event was shifted too far along
    %shift back by median+5;
    ccLag=ccLag(opt.refseis,:);
    Drcm = coralCircShift(Drcm,round(ccLag./Drcm(1).recSampInt));
    evtimes1=datenum([Drcm.recStartTime]');
    [evtimes1 I]=sort(evtimes1);
    Drcm=Drcm(I);
    num=str2num(A(i).name(1:end-4)); %get number of event
    %extract event specific info
    xcorvals=full([Drcm.xcorval]); %xcorvals relative to stack
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
        delete(A(i).name)
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
    save([num2str(num),'.mat'],'Drcm','evtimes1','numev','starttime','endtime','stackAll','recint','recint7','domfreq','f','a','stack7','f7','a7','stack8','f8','a8')
      
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
save('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/CatalogSummary.mat','numevAll','numev7','names',...
    'startdates','enddates','meanRecint','stdRecint','meanRecint7',...
    'stdRecint7','domfreqstack7','signalwidth7','kurtosF7','kurtosT7','evtimes');
end
%% Make set of stacks
%Stack7s=zeros(sizel);
if 1
evtimes=cell(sizel,1);

A=dir('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/*.mat');
for i=1:length(A)
    load(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/IndividualSets/',A(i).name],'stack7','evtimes1');
    num=str2num(A(i).name(1:end-4)); %get number of event
    Stack7s(num)=stack7;
    evtimes{num}=evtimes1;
end
save('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/Stack7s.mat','Stack7s')
save('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/RCMcatalog/evtimes.mat','evtimes')
end
