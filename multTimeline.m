function [allmultTimes]=multTimeline(multiplets,minsize)
% Plot timeline

% Fix MULTTEMPL results to only include multiplet sets with more than
%minsize occurrences and sort in order of first appearance

indSize                = cellfun('size',multiplets,2);
%get rid of sets with minsize or fewer sets
multiplets    = multiplets(indSize>minsize);

%if no events, make blank figure
if isempty(multiplets)
   figure
    return
end

%make it a column vector
if size(multiplets,2)>size(multiplets,1)
    multiplets=multiplets';
end

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
dist=(endTime-startTime)/60;

[junk, I]=sort(dates1);
multiplets=multiplets(I);

clear dates1
clf
hold on
v=1;
allmultTimes=[];

temp=[multiplets{:}];
minim=min([temp.xcorval]);
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
   scatter(Multtimes,ys,100,Xcorval','o','filled');
   caxis([minim 1])

   %caxis([0.65 1])
   text(endTime,v+0.2,num2str(length(Multtimes)),'FontSize',14)
   %multiplet name YYMMDDHHSS of first occurrence
   
   text(startTime+dist,v+0.2,datestr(min(Multtimes),'yymm-ddHHMMSS'));%startTime-dist,v,num2str(i),'FontSize',14)
   
  
   v=v+1;
   allmultTimes=[allmultTimes Multtimes];
   clear Multtimes
    end
end
A=get(gca,'Position');
set(gca,'Position',A-[A(1)/1.5 0 0 0])

for i=1:v-1
     hline(i,':k');
end

datetick('x','mm/dd','keepticks')

hold off
set(gca,'yTick',[]);
%fh = figure(1); % returns the handle to the figure object
set(gcf, 'color', 'white'); % sets the color to white
allmultTimes=sort(allmultTimes);
title(['Timeline of multiplet activity - ',datestr(startTime,'yyyy')],'FontSize',14)
set(gca,'FontSize',14)
ylim([0 i+1])

%put a red line where last update was
if ~isempty(dir('lastupdate.mat'))
    load lastupdate.mat
   vline(lastupdate,'r') 
end

h=colorbar;
xlabel(h,'xcorr','FontSize',14)
B=get(h,'Position');
set(h,'Position',B+[4*B(3) 0 0 0]);
