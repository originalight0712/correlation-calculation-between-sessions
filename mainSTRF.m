% clear
% sessionfolder = 'Elay-230221';
% blockNameexcel = fullfile('D:\Analysis_HC\KilosortPipeline\preprocess4ks@Monty','UsefulBlockNamesEverySession.xlsx');
% [~,~,alldata]=xlsread(blockNameexcel);
% sessRowNum = find(cell2mat(cellfun(@(x) strcmp(x,sessionfolder),alldata(:,1),'UniformOutput', false)));
% blockNames = alldata(sessRowNum,2:end);
% blockNames(cell2mat(cellfun(@(x) any(isnan(x)),blockNames,'UniformOutput', false)))=[];
% 
% % get clusterwise spike in each blocks
% phydataPath = fullfile('D:\Analysis_HC\KilosortPipeline\KS_processedData\',sessionfolder,'KS3_Outp');
% toolboxpath = 'D:\Analysis_HC\toolboxdownload';
% spikeTimeinSampstruct = getSpikeTimeStruct(phydataPath,blockNames,toolboxpath);
% % get strf
% dmrBlock = blockNames(contains(blockNames,'dmr'));
% sprfilepath = ['D:\Analysis_HC\STRF_files\','DNR_Cortex_96k5min_4_50.spr'];
% tankpath = 'F:\Data\';
% savepath = 'D:\Analysis_HC\STRFfig';

clear
sessionfolder = 'Elay-230620';
Ptarray = [];
Pfarray = [];
namearray = [];
Fmarray = [];
RDarray = [];
RTFarray = [];
cc_array = [];
hisarray = [];
blockNameexcel = fullfile('C:\Users\origi\Desktop\summer\strfcreate\strfDemo\strfDemo','UsefulBlockNamesEverySession.xlsx');
[~,~,alldata]=xlsread(blockNameexcel);
sessRowNum = find(cell2mat(cellfun(@(x) strcmp(x,sessionfolder),alldata(:,1),'UniformOutput', false)));
blockNames = alldata(sessRowNum,2:end);
blockNames(cell2mat(cellfun(@(x) any(isnan(x)),blockNames,'UniformOutput', false)))=[];
    
    % get clusterwise spike in each blocks
phydataPath = fullfile('C:\Users\origi\Desktop\summer\strfcreate\strfDemo\strfDemo\data4KS',sessionfolder,'KS3_Outp');
toolboxpath = 'C:\Users\origi\Desktop\summer\strfcreate\strfDemo\strfDemo\toolboxdownload';
spikeTimeinSampstruct = getSpikeTimeStruct(phydataPath,blockNames,toolboxpath);
    % get strf
dmrBlock = blockNames(contains(blockNames,'dmr'));
sprfilepath = fullfile(toolboxpath,'DNR_Cortex_96k5min_4_50.spr');
tankpath = 'C:\Users\origi\Desktop\summer\strfcreate\strfDemo\strfDemo\SynapseTank\Data\';
savepath = 'C:\Users\origi\Desktop\summer\strfcreate\strfDemo\strfDemo\STRFfig1\';
a = 0;
for bb = 1:length(dmrBlock)
    spike_dmr_temp = spikeTimeinSampstruct.(dmrBlock{bb});
    [his,cc,taxis,faxis,STRF1,Fm,RD,Fm1,RD1,name,num,RTFsaving] = getSTRF(sprfilepath,tankpath,sessionfolder,dmrBlock{bb},spike_dmr_temp,fullfile(toolboxpath,'matlab'), savepath);


    Ptarray = [Ptarray,taxis];
    Pfarray = [Pfarray,faxis];
    a = a + num;
    hisarray = [hisarray, his];
    namearray = [namearray,name];
    Fmarray = [Fmarray Fm];
    RDarray = [RDarray RD];
    RTFarray = [RTFarray RTFsaving];
    cc_array = [cc_array cc];

    size_cc = size(cc_array,2);
    depth_array = zeros(1,size_cc/2);
    size1 = size_cc/2;
    for i=1:size1
        depth_array(1,i) = cc_array(1,i*2);
    end
    
end


%% using 3 as basis for the correlation histogram
%ar = [1,2,3,4,5,6,7,8,9,13,15,19,20,21,23,25];% 620
%ar =[1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,18,19,21,23,24,25,26,27,28,30,31];%616
%ar = [1,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25];%613
%ar = [1,3,4,5,6,9,10,11,13,14,15,16,17,18,19,22];%908
%ar = [2,3,5,6,7,9,10,11,12,13,14,15,17];%927
%ar = [1,4,6,7,8,9,10,12,14,16,20,21];%1128
%ar = [1,2,3,4,5,8,10,11];%0922
%ar = [1,3,7,12,15,17,20,26,27,33];%420
%ar = [1,2,3,4,5,6,11,13,14,15,16,17,18,19,20,21,22,23];%705
%ar = [2,4,5,6,10,12,13,14];%711
%ar = [2,3,4,5,6,9,10,11,12,13,14,15];%717
%ar = [1,2,3,8,9,10,11,12,13];
%ar = [2,3,4,5,6,8,9,10,11];%915
%% 
ar = [1,2,3,4,5,6,7,8,9,13,15,19,20,21,23,25];% 620
selectedarray = depth_array(ar);
[B,I] = sort(selectedarray);

size1 = size(ar,2);
Ptmatrix = zeros(size1,size1);
n = size(Ptarray,2);
n1 = n/a;
Pttotal = zeros(n1,size1);

for i=1:size1
    num1 = I(i);
    num2 = ar(num1);
    
    Pttotal(:,i) = Ptarray(num2*n1-n1+1:num2*n1);
 

end

for i=1:size1
    num1 = I(i);
    tmp1 = i;
    for j=1:size1
        num2 = I(j);
        tmp2 = i;
        r = corrcoef(Pttotal(:,tmp1),Pttotal(:,tmp2));
        Ptmatrix(i,j) = r(1,2);
    end
end

for i=1:size1
    for j=1:size1
        r = corrcoef(Pttotal(:,i),Pttotal(:,j));
        Ptmatrix(i,j) = r(1,2);
    end
end

%% 

num_upper = 0;
for i=1:size1
    if(B(i)<13)
        num_upper = num_upper + 1;
    end
end

pt1 = Ptmatrix(1:num_upper,1:num_upper);
pt2 = Ptmatrix(num_upper:size1,num_upper:size1);



figure;
imagesc(Ptmatrix);
axis square
xticks(1:a);
xticklabels(B);
yticks(1:a);
yticklabels(B);
colormap jet;
caxis([0,1]);
colorbar;
title('correlation of time distribution');
annotation('textbox',[0 0.45 0.2 0.4],'string',{['upper average corr ', num2str(round(mean(mean(pt1)),2)) ],...
        ['lower average corr',num2str(round(mean(mean(pt2)),2))],...
        ['average corr ',num2str(round(mean(mean(Ptmatrix)),2))],...
        ['standard deviation',num2str(round(std(Ptmatrix,0,'all'),2))],...
        ['standard deviation of upper',num2str(round(std(pt1,0,'all'),2))],...
        ['standard deviation of lower',num2str(round(std(pt2,0,'all'),2))],...
        },'EdgeColor','none');
saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of time distribution.png']))



%% for peakdelay
Peakdelaymatrix = zeros(size1,size1);
k = size(Fmarray,2);
n2 = k/a;
Peakdelaytotal = zeros(n2,size1);
for i=1:size1
    num1 = I(i);
    num2 = ar(num1);
    Peakdelaytotal(:,i) = Fmarray(num2*n2-n2+1:num2*n2);
        
end
max_value = max(Peakdelaytotal);
figure;
plot(Peakdelaytotal);

xticklabels(B);
xticks(1:size1);
xlabel('channel of the neuron')
ylabel('peak frequency/hz')
saveas(gcf,fullfile(savepath,[sessionfolder,'result of best frequency.png']))

%% 
num_upper = 0;
for i=1:size1
    if(B(i)<13)
        num_upper = num_upper + 1;
    end
end

pf1 = mean(mean(Peakdelaymatrix(1:num_upper,1:num_upper)));
pf2 = mean(mean(Peakdelaymatrix(num_upper:size1,num_upper:size1)));


figure;
imagesc(Peakdelaymatrix);
axis square
caxis([0,1]);
xticks(1:a);
xticklabels(B);
yticks(1:a);
yticklabels(B);
title('correlation of frequency distribution')
colormap jet;
colorbar;
annotation('textbox',[0 0.45 0.2 0.4],'string',{['upper average corr ', num2str(round(mean(mean(pf1)),2)) ],...
        ['lower average corr',num2str(round(mean(mean(pf2)),2))],...
        ['average corr ',num2str(round(mean(mean(Peakdelaymatrix)),2))],...
        },'EdgeColor','none');
saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of frequency distribution.png']))

%% for frequency distribution
Pfmatrix = zeros(size1,size1);
k = size(Pfarray,2);
n2 = k/a;
Pftotal = zeros(n2,size1);
for i=1:size1
    num1 = I(i);
    num2 = ar(num1);
    Pftotal(:,i) = Pfarray(num2*n2-n2+1:num2*n2);
    
    
end
for i=1:size1
    num1 = ar(i);
    
    for j=1:size1
        num2 = ar(j);
        r = corrcoef(Pftotal(:,i),Pftotal(:,j));
        Pfmatrix(i,j) = r(1,2);
    end
end
%% 
num_upper = 0;
for i=1:size1
    if(B(i)<13)
        num_upper = num_upper + 1;
    end
end

pf1 = Pfmatrix(1:num_upper,1:num_upper);
pf2 = Pfmatrix(num_upper:size1,num_upper:size1);


figure;
imagesc(Pfmatrix);
axis square
caxis([0,1]);
xticks(1:a);
xticklabels(B);
yticks(1:a);
yticklabels(B);
title('correlation of frequency distribution')
colormap jet;
colorbar;
annotation('textbox',[0 0.45 0.2 0.4],'string',{['upper average corr ', num2str(round(mean(mean(pf1)),2)) ],...
        ['lower average corr',num2str(round(mean(mean(pf2)),2))],...
        ['average corr ',num2str(round(mean(mean(Pfmatrix)),2))],...
        ['standard deviation of upper',num2str(round(std(pf1,0,'all'),2))],...
        ['standard deviation of lower',num2str(round(std(pf2,0,'all'),2))],...
        },'EdgeColor','none');
saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of frequency distribution.png']))

%% 
num_y = size(RTFarray,2);
num_x = size(RTFarray,1);
num_RTF = num_y/a;
average_RTF = zeros(num_x,num_RTF);

for i=1:a
    average_RTF = RTFarray(:,(i-1)*num_RTF+1:i*num_RTF) + average_RTF;
end
final_RTF = average_RTF/a;
figure;
imagesc(Fm1,RD1,final_RTF),colormap jet
hold on
scatter(Fmarray,-RDarray,[],'black');
hold off;
saveas(gcf,fullfile(savepath,[sessionfolder,'average of RTF.png']))
%% 
RTFtotal = zeros(num_x,num_RTF,size1);
RTFmatrix = zeros(size1,size1);
for i=1:size1
    num1 = I(i);
    num2 = ar(num1);

    for j=1:size1
        num3 = I(j);
        num4 = ar(num3);
        

        r = corrcoef(RTFarray(:,(num2-1)*num_RTF+1:num2*num_RTF),RTFarray(:,(num4-1)*num_RTF+1:(num4)*num_RTF));
        RTFmatrix(i,j) = r(1,2);
    end
end

num_upper = 0;
for i=1:size1
    if(B(i)<13)
        num_upper = num_upper + 1;
    end
end

RTF1 = mean(mean(RTFmatrix(1:num_upper,1:num_upper)));
RTF2 = mean(mean(RTFmatrix(num_upper:size1,num_upper:size1)));
mean(mean(RTFmatrix))

%% 
figure;
imagesc(RTFmatrix);
axis square
xticks(1:a);
xticklabels(B);
yticks(1:a);
yticklabels(B);
colormap jet;
caxis([0,1]);
colorbar;
title('correlation of RTF matrix');
annotation('textbox',[0 0.45 0.2 0.4],'string',{['upper average corr ', num2str(round(mean(mean(RTF1)),2)) ],...
        ['lower average corr',num2str(round(mean(mean(RTF2)),2))],...
        ['average corr ',num2str(round(mean(mean(RTFmatrix)),2))],...
        },'EdgeColor','none');
saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of RTF matrix.png']))

%% histogram
selectedarray = depth_array(ar);
[B,I] = sort(selectedarray);

size1 = size(ar,2);

num_y = size(hisarray,2);
num_x = size(hisarray,1);
num_his = num_y/a;
average_his = zeros(num_x,num_his);
for i=1:a
    average_his = hisarray(:,(i-1)*num_his+1:i*num_his) + average_his;
end
final_his = average_his/a;
figure;
imagesc(final_his),colormap jet
hold on

saveas(gcf,fullfile(savepath,[sessionfolder,'average of histogram.png']))

%% 
histotal = zeros(num_x,num_his,size1);
hismatrix = zeros(size1,size1);
for i=1:size1
    num1 = I(i);
    num2 = ar(num1);

    for j=1:size1
        num3 = I(j);
        num4 = ar(num3);
        

        r = corrcoef(hisarray(:,(num2-1)*num_his+1:num2*num_his),hisarray(:,(num4-1)*num_his+1:(num4)*num_his));
        hismatrix(i,j) = r(1,2);
    end
end

num_upper = 0;
for i=1:size1
    if(B(i)<13)
        num_upper = num_upper + 1;
    end
end

his1 = mean(mean(hismatrix(1:num_upper,1:num_upper)));
his2 = mean(mean(hismatrix(num_upper:size1,num_upper:size1)));
mean(mean(hismatrix))

%% 
figure;
imagesc(hismatrix);
axis square
xticks(1:a);
xticklabels(B);
yticks(1:a);
yticklabels(B);
colormap jet;
caxis([0,1]);
colorbar;
title('correlation of RTF matrix');
annotation('textbox',[0 0.45 0.2 0.4],'string',{['upper average corr ', num2str(round(mean(mean(his1)),2)) ],...
        ['lower average corr',num2str(round(mean(mean(his2)),2))],...
        ['average corr ',num2str(round(mean(mean(hismatrix)),2))],...
        },'EdgeColor','none');
saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of Histogram matrix.png']))
%% useless
%[B,index] = sort(Ptmatrix(1,:),'descend');
name1 = [];
name2 = [];

for i=1:size1
    num1 = ar(i);
    name1 = [name1 num1];
    name2 = [name2 Fmarray(num1)];

    for j=1:size1

        r = corrcoef(Pttotal(:,i),Pttotal(:,j));
        Ptmatrix(i,j) = r(1,2);
    end
end

figure;
imagesc(Ptmatrix);
axis square
xticks(1:a);
xticklabels(name1);
xlabel('depth')
yticks(1:a);
yticklabels(name2);
ylabel('PeakDelay')
colormap jet;
colorbar;
title('correlation of time distribution');

saveas(gcf,fullfile(savepath,[sessionfolder,'correlation of time distribution.png']))






