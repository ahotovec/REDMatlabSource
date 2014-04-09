function [uniqdetect numstadet ndetect] = combineDetections(datecell,rsec)
% Eliminate redundant event times and mark how many stations and which each
% event was detected on
% allstadt.k@gmail.com, written Fall 2012ish

%INPUTS
%datecell: 1xN cell array, where N is the number of stations/channels, where
%each cell is a 6XM matrix of detection times at each station
%rsec: the maximum time in seconds between event detections at the various
%stations to be considered the same event and for those detections to be
%clustered together

%OUTPUTS
%uniqdetect: 6xN double of unique event times between stations
%numstadet: number of stations that detected each uniqdetect event
%ndetect: a 1xN cell array the same size as datecell saying how many
%stations detected each event time (0 means the detection is redundant)

% %EXAMPLE
% 
% %build datecell matrix for 4 stations
% 
% datecell{1,1}=[2012 9 9 0 0 0; 2012 9 9 0 0 10;2012 9 9 0 0 20];
% datecell{1,2}=[2012 9 9 0 0 1; 2012 9 9 0 0 15;2012 9 9 0 0 25];
% datecell{1,3}=[2012 9 9 0 0 1.5; 2012 9 9 0 0 11;2012 9 9 0 0 26];
% rsec=2;
% [uniqdetect numstadet ndetect] = combineDetections(datecell,rsec);
% 
% %output should be
% shouldbe=[2012 9 9 0 0 0; 2012 9 9 0 0 10;2012 9 9 0 0 25;2012 9 9 0 0 15;2012 9 9 0 0 20];
% numev=[3 2 2 1 1];
% ndetectshould={3 2 1}{0 1 1}{0 0 0};

%Make a ndetect matrix for each station that is the same length as datecell
%ndetect says how many stations had a detection for each event
for n=1:length(datecell)
    ndetect{1,n}=ones(1,size(datecell{1,n},2));
end
    
W   = nchoose(datecell); %returns all the possible combinations of one or more elements of the set datecell
ind = cellfun('size',W,2);

%sort so that the intersections of the most events start first
[junk,I] = sort(ind,'descend');
ind     = ind(I);
W       = W(I);

W   = W(ind>1); %remove iterations that are just one station

for k = 1:size(W,1)
    
    %re-sort D so that the smallest set is intersected first, in order to
    %reduce the size of the input matrices into timediffMat and speed the
    %process
    D       = W{k};
    L       = cellfun('size',D,2); %size of inside of W (how many events at each)
    [junk,ind] = sort(L,'ascend');
    D       = D(ind);
    
    %find unique dates that have level detections, and save the ones
    %that should be set to zero because their detection within rsec is already
    %accounted for at another time
    DelT=[];
    level=size(D,2);
    for n = 1:level        
        %DelT=[];        
        if(n==1), Dn = D{n};
        else           
            [junk, I, J] = timediffMat(D{n},Dn,rsec);
            
            %Find the maximum and minimum of the pairs of events that are within rsec between stations and put the maximums
            %in DelT to be deleted and the minimums in Dn to be kept
            ONE=datenum(D{n}(:,I)');
            TWO=datenum(Dn(:,J)');
            DelT=[DelT datevec(max([ONE TWO],[],2))']; %redundant ones to eliminate
            Dn=datevec(min([ONE TWO],[],2))'; %ones to keep reset as Dn, others don't cut it for this level 
        end;  
    end;
    %Change ndetect vector under Dn events to level and DelT events to 0
    for n=1:length(datecell)
        %find events that are in Dn
       [junk,junk,ib]=intersect(datenum(Dn'),datenum(datecell{1,n}'));
       for p=1:length(ib)
           if ndetect{1,n}(ib(p))==1;
               ndetect{1,n}(ib(p))=level;
           end
       end
       [junk,junk,ib]=intersect(datenum(DelT'),datenum(datecell{1,n}'));
       for p=1:length(ib)
           if ndetect{1,n}(ib(p))<=1
               ndetect{1,n}(ib(p))=0;
           end
       end
    end
    
end;

uniqdetect=[datecell{:}];
numstadet=[ndetect{:}];
%get rid of redundant event times
I=numstadet~=0;
numstadet=numstadet(I);
uniqdetect=uniqdetect(:,I);

end


