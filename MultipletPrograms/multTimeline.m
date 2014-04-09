function [allmultTimes]=multTimeline(multiplets,minsize)
% Plot timeline

% Fix MULTTEMPL results to only include multiplet sets with more than
%minsize occurrences and sort in order of first appearance

indSize                = cellfun('size',multiplets,2);
%get rid of sets with 10 or fewer sets
multiplets    = multiplets(indSize>minsize);

%sort by the earliest multiplet date
clear dates1
for i=1:size(multiplets,1)
    dates1(i)=1e8;
    dates2(i)=0;

    for j=1:size(multiplets{i,1},2)       
        if datenum(multiplets{i,1}(1,j).recStartTime')<dates1(i)
         dates1(i)=datenum(multiplets{i,1}(1,j).recStartTime');
        end
        if datenum(multiplets{i,1}(1,j).recStartTime')>dates2(i)
            dates2(i)=datenum(multiplets{i,1}(1,j).recStartTime');
        end
    end
    
end

startTime=min(dates1);
endTime=max(dates2);
dist=(endTime-startTime)/30;

[~, I]=sort(dates1);
multiplets=multiplets(I);

%figure(2)
clear dates1
clf
hold on
v=1;
allmultTimes=[];


for i=1:size(multiplets,1)
    if size(multiplets{i,1},2)>minsize
    for j=1:size(multiplets{i,1},2)
   Multtimes(j)=datenum(multiplets{i,1}(1,j).recStartTime');
   Xcorval(j)=multiplets{i,1}(1,j).xcorval;
    end
   [Multtimes I]=sort(Multtimes); 
   Xcorval=Xcorval(I);
   ys=ones(size(Multtimes));
   ys=ys*v;
   scatter(Multtimes,ys,100,Xcorval,'o','filled');
   %caxis([0.65 1])
   text(endTime+dist,v,num2str(length(Multtimes)));
   text(startTime-dist,v,num2str(i));
   
  
   v=v+1;
   allmultTimes=[allmultTimes Multtimes];
   clear Multtimes
    end
end
for i=1:v-1
     hline(i,':k');
end

datetick('x','mm/dd','keepticks')


hold off
set(gca,'yTick',[]);

allmultTimes=sort(allmultTimes);
title(['timeline of multiplet activity - ',datestr(startTime,'yyyy')])