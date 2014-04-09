multiplets=multtempl.multiplets;
minsize=10;
sum1=1;
    
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

[~, I]=sort(dates1);
multiplets=multiplets(I);

for i=1:length(multiplets)

    if length(multiplets{i})>10
        [cc,ind]=raCrossCorr(multiplets{i});
        A=coralCircShift(multiplets{i},ind(:,end)); 
        figure(1)
        clf;coralPlot(A(1,1:10));
        title([A(1,1).staCode,' - Multiplet ',num2str(sum1),' - ',num2str(length(multiplets{i})),' events'])
        xlabel('seconds')
        saveas(gcf,['multiplet',num2str(sum1),A(1,1).staCode,'.png'])
        
        figure(2)
        clf;coralPlot(A(2,1:10));
        title([A(2,1).staCode,' Multiplet ',num2str(i),' - ',num2str(length(multiplets{i})),' events'])
        pause
        saveas(gcf,['multiplet',num2str(sum1),A(2,1).staCode,'.png'])
        sum1=sum1+1;
    end
    
     
    
end

%%
for i=1:length(multtempl.template)

   clf;coralPlot(multtempl.template{i}) 
   title(['Multiplet ',num2str(i),' - ',num2str(length(multtempl.multiplets{i})),' events'])

         
   pause 
    
end 

%%

 A=multtempl.multiplets{8,1};
% for i=1:100:4000;
%   clf;coralPlot(A(i:i+40));shg
%   pause
%     
% end