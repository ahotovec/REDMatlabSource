function monthvec=monthvector(startYear,startMonth,endYear,endMonth)
%function mon=monthvec(startYear,startMonth,endYear,endMonth)
%make a vector of the start of each month in datenum format

monthvec=[];
for j=startYear:endYear
    if j==startYear && j~=endYear
        addl=[];
        for k=startMonth:12;
            addl=[addl datenum([j k 1 0 0 0])];
        end 
    elseif j==endYear && j~=startYear
        addl=[];
        for k=1:endMonth
            addl=[addl datenum([j k 1 0 0 0])];
        end
    elseif j==endYear && j==startYear
        addl=[];
        for k=startMonth:endMonth
            addl=[addl datenum([j k 1 0 0 0])];
        end
    else
        addl=[];
        for k=1:12
            addl=[addl datenum([j k 1 0 0 0])];
        end
    end
    monthvec=[monthvec addl];
end
