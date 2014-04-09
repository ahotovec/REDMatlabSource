function [D]=coralRemoveDuplicates(D,timebuffer)
%removes any events that appear in more than one set in the cell array of
%structures, keeps the event in the set where it has a higher xcorr value

%INPUTS
%D: cell array of coral structures where each row contains a set of
%repeaters, these must contain the xcorval field
%timebuffer: the time, in seconds, between start time of events to ensure they
%are distinct events (5 sec is a good start)

%OUTPUTS
%D: same as input D but with repeats removed


for k=1:length(D)

    if ~isempty(D{k})
        for n=k+1:length(D) %compare times against each other set
            if ~isempty(D{n}) && ~isempty(D{k})
                times1=datenum([D{k}.recStartTime]');
                times2=datenum([D{n}.recStartTime]');
                A=repmat(times1,1,length(times2));
                B=repmat(times2',length(times1),1);
                differ=abs(A-B)*86400; %find differences in seconds
                rm1=[];rm2=[];
                if min(differ(:))<timebuffer
                    [row col]=find(differ<timebuffer); %find start times within 5 seconds of each other
                    for b=1:length(row) %take out the one with the lower xcorvalue
                        if D{k}(row(b)).xcorval<D{n}(col(b)).xcorval
                            rm1=[rm1 row(b)];
                        else
                            rm2=[rm2 col(b)];
                        end
                    end
                    rm1=unique(rm1);rm2=unique(rm2);
                    D{k}(rm1)=[];
                    D{n}(rm2)=[];
                end
            end
        end
    end
end