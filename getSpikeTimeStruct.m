function spikeTimeinSampstruct=getSpikeTimeStruct(phydataPath,blockNames,phytoolboxpath)
% get spike time in sample, individual block data saved in separate fields, each
% field includes spike times for every clsuter in separate subfields. 
%Spike time in each block are offset so that spike time in each block align to the start of the block. 

% clear
% phytoolboxpath = 'D:\Analysis_HC\toolboxdownload';
% sessiondataPath = 'D:\Analysis_HC\KilosortPipeline\KS_processedData\Elay-230221';
% phydataPath = fullfile(sessiondataPath,'KS3_Outp');
% 
% blockNameexcel = fullfile('D:\Analysis_HC\KilosortPipeline\preprocess4ksCode@Monty','UsefulBlockNamesEverySession.xlsx');
% [~,~,alldata]=xlsread(blockNameexcel);
% blockNames = alldata(1,2:end);


SF_neuro = 24414.0625;
disp('WARNING: assuming spike times are sampled at 24414.0625Hz in the raw.sev files!')
%set path
addpath(genpath(fullfile(phytoolboxpath,'npy-matlab')));
addpath(genpath(fullfile(phytoolboxpath,'spikes')));

% read data
spTime = readNPY(fullfile(phydataPath,'spike_times.npy'));
clID = readNPY(fullfile(phydataPath,'spike_clusters.npy'));
%get cluster info
Clinfo = getClusterInfo_KS3(phydataPath);

% get sample points in each block: ws dmr task1... ws dmr
fnOut_c = fullfile(fileparts(phydataPath),'blockDataPoints.dat');
fidOut_c = fopen(fnOut_c,'r');
blocksamps = fread(fidOut_c,inf,'int64');
fclose('all');

% get clusterwise spike time for each block
spikeTimeinSampstruct = cell2struct(repmat({[]},numel(blockNames),1),blockNames,1);
spikeTimeinSampstruct.fs = SF_neuro;
spikeTimeinSampstruct.Clsinfo = Clinfo;
sampStart=1;
for bb=1:length(blockNames)    
    sampEnd = sum(blocksamps(1:bb));
    indwin = find(spTime>=sampStart&spTime<=sampEnd);
    spTime_temp = spTime(indwin)-sampStart+1;% offset spike time in each block to align to the block onset 
    clID_temp = clID(indwin);
    spikeTimeinSampstruct.(blockNames{bb})=clusterwiseSpike(spTime_temp,clID_temp,Clinfo,blocksamps(bb));
    sampStart = 1+sampEnd;
end
end

function spTime_temp_cls=clusterwiseSpike(spTime_temp,clID_temp,Clinfo,blocksamps)
% cls = sort(unique(clID_temp),'ascend');
% ch = Clinfo.ch(ismember(Clinfo.cluster_id,cls));
% manualLabel = Clinfo.group(ismember(Clinfo.cluster_id,cls));
cls = Clinfo.cluster_id;
ch = Clinfo.ch;
manualLabel = Clinfo.group;
for cc=1:length(cls)
    spTime_temp_cls.(['cls',num2str(cls(cc)),'_ch',num2str(ch(cc)),'_',manualLabel{cc}]) = spTime_temp(clID_temp==cls(cc));
end
spTime_temp_cls.blockSamps = blocksamps;
end