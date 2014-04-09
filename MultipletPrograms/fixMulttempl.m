function [multtempl]=fixMulttempl(multtempl,varargin)
% Fix MULTTEMPL results by only keeping those with the minimum number of
% events specificed by minevents, sorting them by the order in which they appear
% and browsing through the remaining to
% delete garbage events
%usage [multtempl]=fixMulttempl(multtempl,minevents,minxcor)
%minevents,minxcor

%delete all events that have less than the specified minimum xcor value
%wiht the stack
if (length(varargin)>=2)
    minxcor=varargin{2};
    for i=1:size(multtempl.multiplets,1)
        I=[multtempl.multiplets{i}.xcorval]>minxcor;
        multtempl.multiplets{i}=multtempl.multiplets{i}(I);
        multtempl.clustIndex{i}      = multtempl.clustIndex{i}(I);
    end
end

%load the desired multtempl version.
if (length(varargin)>=1)
    minevents=varargin{1};
indSize                = cellfun('size',multtempl.multiplets,2);
%get rid of sets with fewer than minevents sets
multtempl.multiplets    = multtempl.multiplets(indSize>=minevents);
multtempl.template      = multtempl.template(indSize>=minevents);
multtempl.clustIndex    = multtempl.clustIndex(indSize>=minevents);
end

% %sort by the earliest multiplet date
% for i=1:size(multtempl.multiplets,1)
%     dates(i)=1e8;
% 
%     for j=1:size(multtempl.multiplets{i,1},2)       
%         if datenum(multtempl.multiplets{i,1}(1,j).recStartTime')<dates(i)
%          dates(i)=datenum(multtempl.multiplets{i,1}(1,j).recStartTime');
%         end
%     end
%     
% end
% 
% [~, I]=sort(dates);
% multtempl.multiplets=multtempl.multiplets(I);
% multtempl.template=multtempl.template(I);
% multtempl.clustIndex=multtempl.clustIndex(I);

%Now plot each one and delete those that are garbage
v=1;
for i=1:size(multtempl.template,1)
    figure(1)
    clf;coralPlot(multtempl.template{i});
    title([num2str(size(multtempl.multiplets{i},2)),' events'])
    figure(2)
    Randevents=round(rand(1,15)*size(multtempl.multiplets{i},2));
    Randevents(Randevents==0)=1; Randevents(Randevents>size(multtempl.multiplets{i},2))=size(multtempl.multiplets{i},2);
    clf;coralPlot(multtempl.multiplets{i}(unique(Randevents)));
    title([num2str(size(multtempl.multiplets{i},2)),' events - ',num2str(mean([multtempl.multiplets{i}.xcorval])),' median xcorr coeff'])
    button=input('press x to delete, any other key to keep\n','s');
    figure(1);figure(2);
    if ~strcmpi(button,'x')
    index(v)=i;
    v=v+1;
    end
end

if exist('index')
%get rid of deleted events
multtempl.multiplets=multtempl.multiplets(index);
multtempl.template=multtempl.template(index);
multtempl.clustIndex=multtempl.clustIndex(index);
else
   multtempl.multiplets=[];
   multtempl.template=[];
   multtempl.clustIndex=[];
end

save('MULTTEMPL_GOOD.mat','multtempl')