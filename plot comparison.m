label = {'0420','0613','0616','0620','0715','0717','0908','1128'};
%upper_corr = [0.88,0.93,0.94,0.88,0.95,0.75,0.71,0.87];
%lower_corr = [0.48,0.94,0.93,0.89,0.75,0.6,0.53,0.66];

upper_corr = [0.68,0.89,0.6,0.92,0.87,0.93,0.86,0.91];
lower_corr = [0.16,0.77,0.59,0.84,0.61,0.57,0.45,0.43];

%upper_corr = [0.47,0.53,0.27,0.75,0.38,0.48,0.56,0.74];
%lower_corr = [0.26,0.37,0.24,0.61,0.28,0.6,0.22,0.51];
figure
for i = 1:8
    
    xticks(0:2)
    plot([1 2],[upper_corr(i) lower_corr(i)]);
    
    set(gca,'Xticklabel',{'average lower_corr','average upper_corr'})
    hold on
    
    
end
xlabel('upper and lower channels')
ylabel('average correlation')
title('average correlation comparison of upper and lower channels')