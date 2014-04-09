for i=1:length(multtempl.multiplets)
    if length(multtempl.multiplets{i})>10
        [cc,ind]=raCrossCorr(multtempl.multiplets{i});
        A=coralCircShift(multtempl.multiplets{i},ind(:,end));
   clf;coralPlot(A(1:10)) 
   title(['Multiplet ',num2str(i),' - ',num2str(length(multtempl.multiplets{i})),' events'])
   saveas(gcf,['multiplet',num2str(i),'.png'])
    end 
        
   pause 
    
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