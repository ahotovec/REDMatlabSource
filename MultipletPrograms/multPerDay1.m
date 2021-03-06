function [tvecT total]=multPerDay1(multiplets,minSize,startTime,endTime)
%calculate the number of multiplets that occur per day
%(multtempl.multiplets is the input)
allmultTimes=[];
datemin=1e7;
for i=1:size(multiplets,1)
    if size(multiplets{i,1},2)>minSize
        for j=1:size(multiplets{i,1},2)
            Multtimes(j)=datenum(multiplets{i,1}(1,j).recStartTime');
            if Multtimes(j)<datemin; datemin=Multtimes(j); end
            
        end
        allmultTimes=[allmultTimes Multtimes];
        clear Multtimes
    end
end

if numel(allmultTimes)==0
tvecT=startTime:endTime;
total=zeros(size(tvecT));
else

allmultTimes=sort(allmultTimes);

total=zeros(size(startTime:endTime));

   v=1;
    for k=startTime:endTime
       count=allmultTimes(allmultTimes>k & allmultTimes<(k+1));
       count(count>0)=1;
       total(v)=sum(count); 
       v=v+1;
    end

tvecT=startTime:endTime;
% %figure(2)
% plot(tvecT,total)
% datetick('x','mm/dd-HH:MM','keepticks')
% title('Repeating earthquakes per hour')
% hold off 
end
%saveas(gcf,['MULTPERHR',num2str(dayvec(1)),'-',num2str(dayvec(end))],'png')
%save('multperhr.mat','tvecT','total')