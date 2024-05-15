

function C = getClusterInfo_KS3(pathway)
% get cluster information and replace KS group with manual modified group 
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
% pathway = 'F:\Huaizhen\KilosortPipeline\rawContiData\230221\KS3_Outp\';

filename_Autoinfo = fullfile(pathway,'cluster_info.tsv');
C = readtable(filename_Autoinfo, 'Delimiter', '\t','filetype','text');

% filename_Manualinfo = fullfile(pathway,'cluster_group.tsv');
% manC = readtable(filename_Manualinfo, 'Delimiter', '\t','filetype','text');
% C.group = manC.group;


% numerical label for manual group
isUn = cellfun(@(x)strcmp(x,'unsorted'),C.group);
isNoise = cellfun(@(x)strcmp(x,'noise'),C.group);
isMUA = cellfun(@(x)strcmp(x,'mua'),C.group);
isGood = cellfun(@(x)strcmp(x,'good'),C.group);
isDrift = cellfun(@(x)strcmp(x,'drift'),C.group);

cgs = zeros(size(C.group));

cgs(isGood) = 1;
cgs(isMUA) = 2;
cgs(isNoise) = 3;
cgs(isDrift) = 4;
cgs(isUn) = 5;
% add numerical label to cluster info as a new column
C = addvars(C,cgs,'NewVariableNames','numGroup');

