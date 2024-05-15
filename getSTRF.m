%% Load trigger file and spike ties (in seconds) for each cluster
function [hisarray,cc_array,Ptarray,Pfarray,STRF1sarray,Fm1,RD1,Fm3,RD3,namearray,num,RTFsaving] = getSTRF(sprfilepath,tankpath,sessionfolder,Block,spikeTimeStruct,kecktoolboxpath, savepath)
addpath(genpath(kecktoolboxpath));
%% create folder to save strf images and data
if ~exist(fullfile(savepath,sessionfolder),'dir')
    mkdir(fullfile(savepath,sessionfolder))
end
savepath = fullfile(savepath,sessionfolder);

%% Some constants
clusterfieldnames = fieldnames(spikeTimeStruct);
%% get DMR stim Trig times 
%% only works in recording computer because of activecontrol 
% [Data] = readtank_mwa_input_tb(fullfile(tankpath,sessionfolder),Block,1,'local','eSPK');
%% try to read from raw dmr file, low sampling rate, lead to error when run strf code
% dat_temp = TDTbin2mat(fullfile(tankpath,sessionfolder,Block),'STORE','DMRg','TYPE', {'streams'}).streams.DMRg; 
% triggVec = find(dat_temp.data==1);
% Data.Trig = [triggVec(1),triggVec(find(diff(triggVec)>1)+1)]/dat_temp.fs;
% Data.Fs = dat_temp.fs;
%% save the trigger when running channel strf 
load([fullfile(tankpath,sessionfolder,Block),filesep,'dmr_Trigger.mat'])

TrigTimes=round(Data.Fs*Data.Trig);
%          2 presentations
[TrigA,TrigB]=trigfixstrf2(TrigTimes,1800,899); % test
fs=Data.Fs;
disp(['DMR sampling frequency: ', num2str(fs),'Hz'])
Trig = Data.Trig;
Ptarray=[];
Pfarray=[];
STRF1sarray=[];
namearray = [];
Fm1 = [];
RD1 = [];
cc_array = [];
hisarray = [];
num = 0;
RTFsaving = [];
for cc=1:numel(fieldnames(spikeTimeStruct))
%for cc=1:2
    %% Load spike times per cluster
    % Assuming st_clu cell has
    % Take the spike times for the 1st cluster.
    spet = double(spikeTimeStruct.(clusterfieldnames{cc})');
    cc_temp = clusterfieldnames{cc};
    cc_temp(strfind(cc_temp,'_'))='-';
    
    if length(spet)>600
        %% Calculate STRFs for each cluster and trigger block.
        % Parameters
        T1 = 0;
        T2 = 0.15;
        SPL = 75;
        MdB = 30;
        ModType = 'dB';
        Sound = 'MR';
        NBlocks = 1700;
        UF = 2;
        sprtype='float';        
        num_str = regexp(cc_temp,'\d*\.?\d*','match');
        a = size(num_str,2);
        if a>=1            
            num_ch = str2double(num_str(a));
            cc_array = [cc_array, num_ch];
        end
        % Define figure handler for STRF plot to avoid annoyance.
        %h_strf = figure;
        % Calculate STRFs
        % For block A
        [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(sprfilepath,T1,T2,spet,TrigA,fs,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        % For block B
        [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(sprfilepath,T1,T2,spet,TrigB,fs,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        %Average over Triger A and B
        STRF1 = (STRF1A+STRF1B)/2;
        STRF2 = (STRF2A+STRF2B)/2;
        No1=No1A+No1B;
        Wo1=(Wo1A+Wo1B)/2;
        No2=No2A+No2B;
        Wo2=(Wo2A+Wo2B)/2;
        threshold = max(max(STRF1)) * 0.15;
        i_exc = find(STRF1<=threshold); % excitatory part
        i_inh = find(STRF1>=threshold); % inhibitory part
        STRF1e = STRF1; STRF1i = STRF1;
        STRF1e(i_exc) = threshold; STRF1i(i_inh) = threshold;
        %No1 = round((No1A_Ch1 + No1B_Ch1)/2);
        %Wo1 = round((Wo1A_Ch1 + Wo1B_Ch1)/2);
        STRFData(cc) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN,'STRF1e',STRF1e,'STRF1i',STRF1i);
%         RF1P(cc) = strfparam2_hc(STRFData(cc).taxis,STRFData(cc).faxis,STRFData(cc).STRF1e,STRFData(cc).Wo1,STRFData(cc).PP,'MR',500,4,0.5,0.05,0.1,'y',cc_temp);
%add the line
        [STRF1s,Tresh]=wstrfstat(STRFData(cc).STRF1,0.001,STRFData(cc).No1,STRFData(cc).Wo1,STRFData(cc).PP,30,'dB','MR','dB')

        [Fm,RD,RTF,TFParam,RF1P(cc)] = strfparam2_hc(STRFData(cc).taxis,STRFData(cc).faxis,STRF1s,STRFData(cc).Wo1,STRFData(cc).PP,'MR',500,4,0.5,0.05,0.1,'n',cc_temp);
       
        [STRF1sig,Tresh1] = wstrfstat(STRFData(cc).STRF1,0.001,STRFData(cc).No1,STRFData(cc).Wo1,STRFData(cc).PP,30,'dB','MR','dB');
        
        % add the line
       
%         figure('visible','off'); 
%         taxis=(taxis)*1e3;
%         % faxis=(faxis)*1e3;
%         pcolor(taxis,log2(faxis/faxis(1)),STRF1);
%     %         pcolor(taxis,log2(faxis/faxis(1)),STRF1A); % test
%         colormap jet;set(gca,'YDir','normal'); shading flat;
%          %Use this to get BF excitation
%         best_freq = faxis(1) * 2^RF1P(cc).PeakBFP;
%        
        

        %[RDHist1,FMHist1,RDHist2,FMHist2,Time1,Time2]=rtfhist(sprfilepath,spet,TrigA,fs);
        
        
        Ptarray = [Ptarray,RF1P(cc).Pt];
        Pfarray = [Pfarray,RF1P(cc).Pf];
        STRF1sarray = [STRF1sarray,STRF1s];
        [BF,his,RTF1,Fm3,RD3] = plotRTFSTRF(taxis,faxis,Fm,RD,RTF,TFParam,STRF1s,STRF1sig,RF1P(cc),cc_temp,sprfilepath,spet,TrigB,fs)
        C = cellstr(cc_temp);
        namearray = [namearray,C];
        cc_array = [cc_array, num_ch];
        saveas(gcf,fullfile(savepath,[sessionfolder,'_',Block,'_STRF_',cc_temp,'.png']))
                % 创建一个包含数字字段的结构体
        num = num + 1;
        Fm2 = round(RF1P(cc).BestFm(1),2);
        RD2 = round(RF1P(cc).PeakEnvDelay,2);
        Fm1 = [Fm1 Fm2];
        RD1 = [RD1 RD2];
        RTFsaving = [RTFsaving RTF1];
        hisarray = [hisarray,his];
        
        
        
        
%         save(fullfile(savepath,[sessionfolder,'_',Block,'_STRF_',cc_temp,'.mat']),'STRFData')
        %% Significant STRFs
        %p = 0.05;
        %SModType = 'dB';
        %[STRF1s,Tresh1]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
        %catch 
        %    disp(['Neuron',num2str(i),'did not work'])
        %end
    end
end
end
%% 


function [BestFm,N,RTF1,Fm1,RD1] = plotRTFSTRF(taxis,faxis,Fm,RD,RTF,TFParam,STRF1,STRF1sig,RFParam,chStr,sprfilepath,spet,TrigA,fs)
    X=log2(faxis/faxis(1));             %Octabes

    figure('Position',[10 10 900 600],'visible','off')
    p=[ 0.2915-0.08    0.2872-0.08    0.6135    0.6378];
    subplot('Position',p)
    subplot(2,2,1)
    

    
    imagesc(taxis,X,STRF1)    
    colormap jet
    hold on
    %overlay significant contour map
    [C,h] = contourf(taxis,X,STRF1sig,'fill','off');
    set(gca,'YDir','normal')
    %RFParam.PeakBF = fft(RFParam.PeakBF);
    

    hold off
    title(['STRF ',chStr ])
    ylabel('Frequency (octave)')
    xlabel('Delay (s)')
    best_freq = round(faxis(1) * 2^RFParam.PeakBF);            
    annotation('textbox',[0.83 0.45 0.2 0.4],'string',{['BFre ', num2str(best_freq) 'Hz'],...
        ['FBW50%  ',num2str(round(RFParam.HalfEnvBW)),'Hz'],...
        ['Delay ',num2str(round(RFParam.PeakDelay)),'ms'],...
        ['BFm  ',num2str(round(RFParam.BestFm(1),2))],...
        ['BRD  ',num2str(round(RFParam.BestRD(1),2))],...
        },'EdgeColor','none');
    BestFm = RFParam.BestFm(1);
        %subplot(1.2,5,1)


    subplot(2,2,2)
    [RDHist1,FMHist1,RDHist2,FMHist2,Time1,Time2]=rtfhist(sprfilepath,spet,TrigA,fs);
    
    [FMAxis,RDAxis,N]=hist2(-FMHist1,RDHist1,20,20,'n');
    for i=1:20
        for j=1:20
            if(N(i,j)<3)
                N(i,j)=0;
            end
        end
    end
    
    pcolor(FMAxis,RDAxis,N);colormap jet;
    colorbar;
    title('histogram of FM and RD')
    

    
    
    subplot(2,2,3)
    MaxRD=2;
    MaxFm=50;
    [Fm5,RD5,RTF1,RTFVar]=strf2rtf(taxis,faxis,STRF1,MaxFm,MaxRD);
    imagesc(Fm,RD,RTF1),shading flat,colormap jet
    title('RFT figure')
    xlabel('Fm-modulation rate')
    ylabel('ripple density(Cycles/Octave')
    Fm1 = Fm;
    RD1 = RD;
    %STRF2 = fftshift(fft2(STRF1));
    %imagesc(Fm,RD,abs(STRF2)),shading flat,colormap jet
    %title('RFT of forier transform')

    
    %subplot(3,2,4)

    %targetSize = [40, 20];

    %newMatrix = imresize(RTF1, targetSize);
    %correlationValue = corr2(newMatrix, N);
    %text(0.5, 0.5, sprintf('Correlation: %.4f', correlationValue), 'FontSize', 14, 'HorizontalAlignment', 'center');
    %axis off;
    %title('Correlation');


    %subplot(3,2,5)
    %corrematrix = zeros(40,20);    
    %for i=1:40
    %    for j = 1:20
    %        corrematrix(i,j) = (newMatrix(i,j)-N(i,j))/(newMatrix(i,j)+N(i,j));
    %    end
    %end
    %pcolor(FMAxis,RDAxis,corrematrix);colormap jet;
    %colorbar;

    %title('correlation between RTFs and histograms')

    




  

    
    

    




end
%% 
function plotSigStrf(taxis,faxis,STRF1,STRF1sig,RFParam,chStr)
    X=log2(faxis/faxis(1));             %Octabes

    figure('Position',[10 10 900 600],'visible','off')
    p=[ 0.2915-0.08    0.2872-0.08    0.6135    0.6378];
    subplot('Position',p)
    imagesc(taxis,X,STRF1)    
    colormap jet
    hold on
    %overlay significant contour map
    [C,h] = contourf(taxis,X,STRF1sig,'fill','off');
    set(gca,'YDir','normal')
    plot(RFParam.PeakDelay,RFParam.PeakBF,'wx','linewidth',2)
    hold off
    title(['STRF ',chStr ])
    ylabel('Frequency (octave)')
    xlabel('Delay (ms)')
    best_freq = round(faxis(1) * 2^RFParam.PeakBF);            
    annotation('textbox',[0.83 0.45 0.2 0.4],'string',{['BFre ', num2str(best_freq) 'Hz'],...
        ['FBW50%  ',num2str(round(RFParam.HalfEnvBW)),'Hz'],...
        ['Delay ',num2str(round(RFParam.PeakDelay)),'ms'],...
        ['BFm  ',num2str(round(RFParam.BestFm(1),2))],...
        ['BRD  ',num2str(round(RFParam.BestRD(1),2))],...
        },'EdgeColor','none');
    
    
    %subplot(1.2,5,1)
    p=[0.1300-0.08    0.2872-0.08    0.1237    0.6378];
    subplot('Position',p)
    plot(RFParam.Pf/max(RFParam.Pf),X,'k')
    hold on
    plot(0.5,RFParam.X1,'r+')
    plot(0.5,RFParam.X2,'r+')
    plot(0.1,RFParam.XL10,'r.')
    plot(0.1,RFParam.XU10,'r.')
    plot(1,RFParam.PeakEnvBF,'ro')
    axis([0 1.1 0 max(X)])
    set(gca,'XDir','reverse')
    set(gca,'Visible','off')
    hold off
    
    
    %subplot(5,1.2,6)
    p=[0.2915-0.08    0.1100-0.08    0.6135    0.1243];
    subplot('Position',p)
    plot(taxis,RFParam.Pt/max(RFParam.Pt),'k')
    hold on
    plot(RFParam.t1_50,0.5,'r+')
    plot(RFParam.t2_50,0.5,'r+')
    plot(RFParam.t1_10,0.1,'r.')
    plot(RFParam.t2_10,0.1,'r.')
    plot(RFParam.PeakEnvDelay,1,'ro')
    axis([0 max(taxis) 0 1.1])
    set(gca,'YDir','reverse')
    set(gca,'Visible','off')
    hold off
end