%function combineCatalogs
%code to combine RCMcatalog and RCScatalog that were obtained using
%template searching, started writing from combineMultdetect.m program used
%in real-time multiplet detection system, started 10 July 2013

load RCScatalog/evtimes.mat
evT_RCS.dates1=evtimes;
load RCMcatalog/evtimes.mat
evT_RCM.dates1=evtimes;

rsec=10; %maximum time in seconds between detections at the various stations to be considered the same event
perc=0.25; %proportion of detections that must overlap
%%
%make index matrix for each station, find those that match
RCS=1:length(evT_RCS.dates1)';
RCM=1:length(evT_RCM.dates1)';

INDMAT=zeros(length(evT_RCS.dates1),2);
INDMAT(:,1)=RCS;
%%
for i=1:length(evT_RCS.dates1)
    if ~isempty(evT_RCS.dates1{i})
        datecell{1,1}=datevec(evT_RCS.dates1{i})';
        
        for j=1:length(evT_RCM.dates1)
            if RCM(j)==0 || isempty(intersect(evT_RCM.dates1{j}(evT_RCM.dates1{j}>min(evT_RCS.dates1{i})),evT_RCM.dates1{j}(evT_RCM.dates1{j}<max(evT_RCS.dates1{i})))) %only compare if there is some overlap
                continue
            else
                    datecell{1,2}=datevec(evT_RCM.dates1{j})';
                    [uniqdetect numstadet ndetect] = combineDetections(datecell,rsec);
                    I=ndetect{1}==2;
                    %if sum(I)/length(ndetect{1})>0.1 %if more than 10% overlap, they're the same
                    if (length(numstadet(numstadet==2))/length(numstadet))>perc
                        INDMAT(i,2)=RCM(j);
                        RCM(j)=0;
                    end
            end
        end
    end
    
end

%Make a list of the sets that map to each other, left column RCS and right
%column, the corresponding sets in RCM
I=INDMAT(:,2);
I=find(I>0);
Sets=INDMAT(I,:);
STA={'RCS','RCM'};
%save
save('SharedCatalog-RCS-RCM_10s25pc.mat','STA','Sets')



