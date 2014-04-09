%script to extract the template of each multiplet set found (and the number
%of repeats) and then compare all templates in all year with each other

clear all;close all;clc
if 1
A=dir('/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/STARmultsearch');
Temp={};
Num={};
Full={};
%count=1;
for i=4:length(A)
    cd(['/Users/Kate/MATLAB1/Rainier_all/LongTermSearch/STARmultsearch/',A(i).name])
    B=dir('MULTTEMPL_GOOD.mat');
    if ~isempty(B)
        load MULTTEMPL_GOOD.mat
        Temp=[Temp; multtempl.template];
        Full=[Full; multtempl.multiplets];
        if ~isempty(multtempl.multiplets)
        Num=[Num; cellfun(@size,multtempl.multiplets,'UniformOutput',false)];
        end
        %count=count+1;
    end
    
end

for i=1:length(Temp)
   Templates(i)=Temp{i};
   NumEv(i)=Num{i}(2);
end
%extract number of events

cd /Users/Kate/MATLAB1/Rainier_all/LongTermSearch
save('MULTTEMPL_GOOD_TEMPLATES_STAR.mat','Templates','NumEv')
save('MULTTEMPL_GOOD_FULL_STAR.mat','Full','NumEv')
else
    load MULTTEMPL_GOOD_TEMPLATES_STAR.mat
    load MULTTEMPL_GOOD_FULL_STAR.mat
end

%% compare largest sets
if 1

I=find(NumEv>5);
BigT=Templates(I);
NumEvBig=NumEv(I);
BigFull=Full(I);

[cc, t, indx, ccLag, ccMax, ierr] = coralCrossCorr(BigT);
save('StuffSTAR.mat','ccLag','ccMax','BigT','NumEvBig')

end
%% Plot similar waveforms and see how they compare
if 1
takeout=[];
tvec=0:BigT(1).recSampInt:BigT(1).recSampInt*(BigT(1).recNumData-1);
height=0;
for i=1:length(BigT)
    if isempty(intersect(i,takeout))
        I=find(ccMax(i,:)>0.7);
        if length(I)>1
            for j=1:length(I)
                
                plot(tvec,BigT(I(j)).data+height)
                text(0.2,height+0.05,[datestr(BigT(I(j)).recStartTime'),' ',num2str(NumEvBig(I(j))),'ev'],'Color','r')
                hold on
                height=height+0.4;
                
            end
        
        hold off
        axis tight
        
        takeout=[takeout,I];
        
        pause
        end
    end
end
end

%% Group together templates that correlate above 0.8 (that occur within 1
% month of each other)
if 1
    load StuffSTAR.mat
    load MULTTEMPL_GOOD_FULL_STAR.mat
    BigFull=Full;
MasterFull={};%initiate matrix that will hold all waveforms
takeout=[];
count=1;
for i=1:length(BigT)
    if isempty(intersect(i,takeout))
        I=find(ccMax(i,:)>0.8);
        times=datenum([BigT(I).recStartTime]');%get time of first event-ish of each
        %J=find(ccMax(i,:)==max(ccMax(i,:))); %find index of starter event
        timediff=times-datenum(BigT(i).recStartTime');%find time differences between starter event and all other events
        K=I(abs(timediff)<90); %only combine events that are within a few months of each other
        temp=[];
        for j=1:length(K)
            temp=[temp BigFull{K(j)}];
        end
        takeout=[takeout,K];%take out events that have already been combined with others
        MasterFull{count}=temp;
        count=count+1;
    end
end

%get date span of each of them
startdates=zeros(size(MasterFull));
enddates=startdates;
for i=1:length(MasterFull)
   temp1=MasterFull{i};
   temp2=datenum([temp1.recStartTime]');
    startdates(i)=min(temp2);
    enddates(i)=max(temp2);
    
end

save('MASTER_FULL_STAR.mat','MasterFull');

end
%% align all of them 
if 1
    load MASTER_FULL_STAR.mat
TempStack=[];
MasterFull0p7thresh={};
for i=1:length(MasterFull)
    temp=MasterFull{i};
    opt.refseis=1;%only get cross corr relative to first event
    [~, ~, ~, ccLag] = coralCrossCorr(temp, opt);
    %if median shift is greater than 5 seconds, it probably means the first event was shifted too far along
    %shift back by median+5;
    ccLag=ccLag(opt.refseis,:);
    med=median(ccLag);
    if med>5
        ccLag=ccLag-med+5;
    end
    MasterFull{i} = coralCircShift(temp,round(ccLag./MasterFull{i}(1).recSampInt));
    %note the saved version of MASTER_FULL_STAR.mat is the aligned version
    %stack them
    TempStack{i}=coralRaStack(MasterFull{i});
    TempStack{i}=rmfield(TempStack{i},'multiplets');
    %cross correlate each individual event with its stack
    [~, ~, ~, ~, ccMax] = coralCrossCorr([TempStack{i} MasterFull{i}], opt);
    %find individual events that don't correlate above 0.7 and remove them
    ccMax=ccMax(opt.refseis,:);
    I=find(ccMax>0.7)-1; %minus 1 because of addition of TempStack to MasterFull above
    I(I==0)=[];
    MasterFull0p7thresh{i}=MasterFull{i}(I);
    %stack again (assuming there are events that correlate high enough, if
    %not, leaves it empty)
    if ~isempty(I)
    MasterStack{i}=coralRaStack(MasterFull0p7thresh{i});
    MasterStack{i}.multiplets=datenum([MasterFull0p7thresh{i}.recStartTime]');%bug in coralRaStack, replace multiplet times with actual ones
    end
end
%turn into array of structures instead of cell array because cell arrays
%are annoying
if 1
for i=1:length(MasterStack)
    if ~isempty(MasterStack{i})
MasterStack1(i)=MasterStack{i};
    end
end
MasterStack=MasterStack1;

save('MasterStacks_STAR.mat','MasterStack','startdates','enddates')
save('MASTER_FULL_STAR_0.7thresh.mat','MasterFull0p7thresh')
end
end