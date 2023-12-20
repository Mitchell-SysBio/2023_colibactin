%% user variables
curDate = '20230808'; %for saving files at the end
wellInx = 5; %column with well location from metadata
ctrlNum = 2; %number of control samples you have in your metadata/.fcs files
repNum = 3; %number of replicates you have per sample (used for plotting/taking averages)
%% load gating coordinates
load xyCoordinate.mat;
xyCells = xyCoord.Cells; % cells gate
xySingle = xyCoord.Single; % single cells gate
xyCFP = xyCoord.CFP; %CFP gate
xyYFP = xyCoord.YFP; %YFP gate
xymCh = xyCoord.mCherry; % mCherry gate

%% max signal value
maxSigVal = 2^18; %assumed maximum fluorescent signal the machine records

%% read in fcs file names
fcsFiles = dir('*.fcs'); %get all .fcs file names
metadata = readcell('example_metadata.xlsx'); %get metadata descriptions
sampleNum = length(metadata)-ctrlNum; %get the number of samples

%% read files with fca_readfcs
fcsdat = {}; fcshdr = {};
for iFile = 1:size(fcsFiles,1)
    [fcsdat{iFile}, fcshdr{iFile}] = fca_readfcs(fcsFiles(iFile).name); %read in all .fcs files
end
allNames = replace({fcsFiles(:).name},'_',' '); %get fcs file names for all samples

%% sort files into same order as metadata
sampInx = [];
for iWell = 1:size(metadata,1)
    curWell = metadata(iWell,wellInx);
    sampInx(iWell,:) = find(contains(allNames,curWell)); %get index of each file in the order of the metadata
end
orderedFcsFiles = fcsdat(sampInx); %reorder the fcs files to match the metadata
description = join(cellfun(@string,metadata)); %combine the descriptors to one string for each sample
condKey = description([1:repNum:sampleNum]); %get descriptions for averaged samples

%% Apply gating strategy to all fcs files
%extract cell vectors for each parameter for all files
FSCA = {}; FSCH = {}; SSCA = {}; CFP = {}; YFP = {};
for i = 1:length(fcsFiles)
    FSCA{i} = orderedFcsFiles{i}(:,1);
    FSCH{i} = orderedFcsFiles{i}(:,2);
    SSCA{i} = orderedFcsFiles{i}(:,4);
    CFP{i} = orderedFcsFiles{i}(:,8);
    YFP{i} = orderedFcsFiles{i}(:,9);
    mCh{i} = orderedFcsFiles{i}(:,7);
end
%get indexes of cells within each gate
for i = 1:length(fcsFiles)
    tfCells{i} = inpolygon(FSCA{i},SSCA{i},xyCells(:,1),xyCells(:,2)); %get indexes of cells that are inside that gate
    tfSingle{i} = inpolygon(FSCA{i},FSCH{i},xySingle(:,1),xySingle(:,2)); %get indexes of cells inside the gate
    tfYFP{i} = inpolygon(FSCA{i},YFP{i},xyYFP(:,1),xyYFP(:,2)); %get indexes of cells within YFP pos gate
    tfCFP{i} = inpolygon(FSCA{i},CFP{i},xyCFP(:,1),xyCFP(:,2)); %get indexes of cells within CFP pos gate
    tfmCh{i} = inpolygon(FSCA{i},mCh{i},xymCh(:,1),xymCh(:,2)); %get indexes of cells within mch pos gate
end
% get indexes of cells that are within both the cell and single cell gates
tfSingleCells = {};
for i =1:length(fcsFiles)
    tfSingleCells{i} = tfSingle{i} & tfCells{i}; %get indexes of cells within both the cell and single cell gates
    tfSingleCellsCFP{i} = tfSingleCells{i} & tfCFP{i}; %get indexes of single cells in the positive CFP gate
    tfSingleCellsYFP{i} = tfSingleCells{i} & tfYFP{i}; %get indexes of single cells in the positive YFP gate
    tfSingleCellsmCh{i} = tfSingleCells{i} & tfmCh{i}; %get indexes of single cells in the positive mCh gate
end

%% determine number of cells and percent of cells in each gate
%make variable holders
rawYFP = []; rawCFP = []; rawmCh = []; rawSingle = []; histYFP = [];
histCutoff = 10^4.5;
for i = 1:length(fcsFiles)
    rawYFP = [rawYFP, sum(tfSingleCellsYFP{i} & tfSingleCellsCFP{i})]; %sum of YFP positive CFP positive single cells with no normalization
    histYFP = [histYFP, sum(YFP{i}(tfSingleCellsCFP{i})>histCutoff)]; %sum of YFP positive CFP positive single cells by visual cutoff
    rawCFP = [rawCFP, sum(tfSingleCellsCFP{i})]; %sum of CFP positive single cells
    rawmCh = [rawmCh, sum(tfSingleCellsmCh{i})]; %sum of CFP positive single cells
    rawSingle = [rawSingle, sum(tfSingleCells{i})]; %sum of single cells
end

percYFP = (rawYFP./rawCFP)*100; %percent of YFP positive cells
percCFP = (rawCFP./rawSingle)*100; %percent of CFP positive cells
percmCh = (rawmCh./rawSingle)*100; %percent of mcherry positive cells


%% concatenate replicates into a single array
condIdx = [1:3]; % indexes for time points (based on # replicates)

%make variable holders
posTimeYFP = {}; catPosTimeY = {}; catPosY = {}; negTimeYFP = {}; catNegTimeY={}; catNegY = {};
posTimeCFP = {}; catPosTimeC = {}; catPosC = {}; negTimeCFP = {}; catNegTimeC={}; catNegC = {};

for iCond = 1 %number of time points (1:nTime)
    posTimeYFP = YFP(condIdx(iCond,:)); %extract pks+ YFP for time point i
    negTimeYFP = YFP(condIdx(iCond,:)+3); %extract pks- YFP for time point i
    posTimeCFP = FSCA(condIdx(iCond,:)); %extract pks+ CFP for time point i
    negTimeCFP = FSCA(condIdx(iCond,:)+3); %extract pks- CFP for time point i
    for i = 1:length(condIdx(iCond,:))
        catPosTimeY{i} = posTimeYFP{i}(tfSingleCellsCFP{condIdx(iCond,i)}); %extract YFP signal for CFP+ single cells in pks+
        catPosY{iCond} = vertcat(catPosTimeY{:}); %concatenate all YFP values for CFP+ single cells in pks+
        catNegTimeY{i} = negTimeYFP{i}(tfSingleCellsCFP{condIdx(iCond,i)+3}); %extract YFP signal for CFP+ single cells in pks-
        catNegY{iCond} = vertcat(catNegTimeY{:}); %concatenate all YFP values for CFP+ single cells in pks-
        catPosTimeC{i} = posTimeCFP{i}(tfSingleCellsCFP{condIdx(iCond,i)}); %extract CFP signal for CFP+ single cells in pks+
        catPosC{iCond} = vertcat(catPosTimeC{:}); %concatenate all CFP values for CFP+ single cells in pks+
        catNegTimeC{i} = negTimeCFP{i}(tfSingleCellsCFP{condIdx(iCond,i)+3}); %extract CFP signal for CFP+ single cells in pks-
        catNegC{iCond} = vertcat(catNegTimeC{:}); %concatenate all CFP values for CFP+ single cells in pks-
    end
end


%% plot histogram over time with concatenated replicates

posMed = []; negMed = []; %place holders
time = {'12h'}; %string array for time points

figure; hold on;
inx = [1]; %indexes of time points to plot
for i = 1:length(inx)
    subplot(length(inx),1,i); hold on;
    posMed(i) = median(catPosY{inx(i)}); %median YFP signal pks+
    negMed(i) = median(catNegY{inx(i)}); %median YPF signal pks-
    histogram(catPosY{inx(i)},10.^[0:0.025:6],'Normalization','probability','EdgeColor','none') %pks+ YFP histogram
    histogram(catNegY{inx(i)},10.^[0:0.025:6],'Normalization','probability','EdgeColor','none') %pks- YFP histogram
    title(time{inx(i)})
    grid on; box on;
    ylim([0 .1])
    xlim([10^2.5 10^5.5])
    set(gca,'xscale','log')
    xline(posMed(i),'k'); %mark pks+ median
    xline(negMed(i),'r'); %mark pks- median
end
legend({'pks+','pks-'});

%% bar plot % recA positive with cutoff of +3 s.d. from mean of 0h
logPos = log10(catPosY{1}); %convert 0h pks+ to log scale
logNeg = log10(catNegY{1}); %convert 0h pks- to log scale

meanPos = mean(logPos); %mean of 0h pks+
meanNeg = mean(logNeg); %mean of 0h pks-
stdPos = std(logPos); %standard deviation of 0h pks+
stdNeg = std(logNeg); %standard deviation of 0h pks-

histYFP = [];
histCutoff = 10^(mean([meanPos+3*stdPos,meanNeg+3*stdNeg])); %take mean of 0h pks+ and 0h pks- plus 3 s.d.
for i = 1:length(fcsFiles)
    histYFP = [histYFP, sum(YFP{i}(tfSingleCellsCFP{i})>histCutoff)]; %sum of YFP positive CFP positive single cells by cutoff in all replicates and conditions
end

percHistYFP = (histYFP./rawCFP)*100; %percent of YFP positive cells by cutoff for all conditions and replicates

%avg replicates
meanHistYFP = []; stdHistYFP = []; %place holders
for i = [1:repNum:sampleNum]
    meanHistYFP = [meanHistYFP, mean([percHistYFP(i),percHistYFP(i+1),percHistYFP(i+2)])]; %mean of replicates (alternates pks+ and pks-)
    stdHistYFP = [stdHistYFP, std([percHistYFP(i),percHistYFP(i+1),percHistYFP(i+2)])]; %standard deviation of replicates (alternates pks+ and pks-)
end

% plot
time = {'12h'}; %string array for time points
meanHistPlot = [meanHistYFP(1:2:end);meanHistYFP(2:2:end)]'; %extract pks+ mean YFP positive cells in one column and pks- mean YFP positive cells in second column
stdHistPlot = [stdHistYFP(1:2:end); stdHistYFP(2:2:end)]'; %extract pks+ s.d. in column 1, pks- s.d. in column 2

b = bar(1:length(time),meanHistPlot); %bar plot mean YFP positive cells

%add error bars
y = meanHistPlot; %y for each error bar
err = stdHistPlot; %error size for each error bar

% Calculating the width for each bar group to set x for each errorbar
[ngroups,nbars] = size(y); %holders for number of bars
x = nan(nbars, ngroups); %holders for number of bars
for i = 1:nbars
   x(i,:) = b(i).XEndPoints; %get middle of each bar on x axis
end
hold on
errorbar(x', y, err, 'k.'); %plot error bars
hold off

grid on; box on;
set(gca,'xtick',1:length(time),'xticklabel',time);
xlabel('time (hours)');
ylabel('% recA positive reporter cells')
legend({'pks+','pks-'})
title('mean of pks+ and pks- cutoffs')

%% save % positive YFP, CFP+ cells, mCherry+ cells
fileStr = [curDate,'_','recA_results.xlsx'];
mat2save = [percYFP',rawCFP', percCFP', rawmCh', percmCh', rawSingle'];
headerStr = {'condition','% recA positive','#CFP','% CFP','#mCherry','% mCherry','total single cells'};
writecell(headerStr,fileStr,'Range','A1');
writematrix(description,fileStr,'Range','A2');
writematrix(mat2save,fileStr,'Range','B2');
writematrix(condKey,fileStr,'Sheet',2,'Range','A1');
writematrix(meanHistYFP',fileStr,'Sheet',2,'Range','B1');
writematrix(stdHistYFP',fileStr,'Sheet',2,'Range','C1');
