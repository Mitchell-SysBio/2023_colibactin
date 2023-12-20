%% user variables
curDate = '20230808'; %for saving files at the end
wellInx = 5; %column with well location from metadata
ctrlNum = 2; %number of control samples you have in your metadata/.fcs files
repNum = 3; %number of replicates you have per sample (used for plotting/taking averages)

%% gating example indexing
% select representative samples to get positive and negative populations
% for all fluorophores based on metatdata document
unst = 7; %index for unstained ctrl well
recA = 8; %index for recA ctrl well
pos24 = 3; %index for 24h pks+
neg24 = 6; %index for 24h pks-

subset = [recA, pos24,neg24, unst]; %create array of indexes

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

%% Extract vectors for each parameter (indexes can be found in fcshdr.par)
FSCA = {}; FSCH = {}; SSCA = {}; CFP = {}; YFP = {};
for i = 1:length(subset)
    FSCA{i} = orderedFcsFiles{subset(i)}(:,1); %forward scatter area
    FSCH{i} = orderedFcsFiles{subset(i)}(:,2); %forward scatter height
    SSCA{i} = orderedFcsFiles{subset(i)}(:,4); %side scatter area
    CFP{i} = orderedFcsFiles{subset(i)}(:,8); %CFP area
    YFP{i} = orderedFcsFiles{subset(i)}(:,9); %YFP area
    mch{i} = orderedFcsFiles{subset(i)}(:,7); %mCherry area
end
%% Gate on cells from subset samples
myColor = jet(length(subset)); %get as many different colors as there are samples
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i},SSCA{i},'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)); %scatter plot of FSC-A x SSC-A
end
xlabel('FSC-A'); ylabel('SSC-A');
set(gca,'Yscale','log','Xscale','log');
legend(description(subset),'location','eastoutside');
cells = drawpolygon(gca); %draw polygon to gate on cells
xyCells = cells.Position; %get coordinates of the polygon points
close(fh);
tfCells = {};
for i = 1:length(subset)
    tfCells{i} = inpolygon(FSCA{i},SSCA{i},xyCells(:,1),xyCells(:,2)); %get indexes of cells that are inside that gate
end

%% Gate on single cells within the cells gate
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfCells{i}),FSCH{i}(tfCells{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)); %scatter plot of FSC-A x FSC-H
end
xlabel('FSC-A'); ylabel('FSC-H');
set(gca,'Yscale','log','Xscale','log');
legend(description(subset),'location','eastoutside');
singleCell = drawpolygon(gca); %draw polygon to gate on single cells
xySingle = singleCell.Position; %get coordinates
close(fh);
tfSingle = {};
for i = 1:length(subset)
    tfSingle{i} = inpolygon(FSCA{i},FSCH{i},xySingle(:,1),xySingle(:,2)); %get indexes of cells inside the gate
end

%% Get indexes of single cells within cells gate
tfSingleCells = {};
for i =1:length(subset)
    tfSingleCells{i} = tfSingle{i} & tfCells{i}; %get indexes of cells within both the cell and single cell gates
end
%% Gate on CFP positive cells within the single cells gate
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfSingleCells{i}),CFP{i}(tfSingleCells{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)) %scatter plot of CFP in single cells
end
xlabel('FSC-A'); ylabel('CFP-A');
yline(maxSigVal);
legend(description(subset),'location','eastoutside');
xlim([10 100000]);
set(gca,'Yscale','log','Xscale','log');
grid on;
posCFP = drawpolygon(gca); %draw rectangle to gate YFP cells
xyCFP = posCFP.Position; %get coordinates of rectangle
close(fh)
tfCFP = {};
for i = 1:length(subset)
    tfCFP{i} = inpolygon(FSCA{i},CFP{i},xyCFP(:,1),xyCFP(:,2)); %get indexes of cells within CFP pos gate
end
%% Gate on raw YFP positive cells looking at CFP + single + cells
%All CFP+ cells have low YFP, so gate on the brighter cells
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfSingleCells{i}&tfCFP{i}),YFP{i}(tfSingleCells{i}&tfCFP{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)) %scatter plot of YFP signal in single cells
end
yline(maxSigVal);
set(gca,'Yscale','log','Xscale','log');
xlabel('FSC-A'); ylabel('YFP-A');
xlim([10 100000]);
grid on;
legend(description(subset),'location','eastoutside');
posYFP = drawpolygon(gca); %draw rectangle to gate YFP cells
xyYFP = posYFP.Position; %get coordinates of rectangle
close(fh)
tfYFP = {};
for i = 1:length(subset)
    tfYFP{i} = inpolygon(FSCA{i},YFP{i},xyYFP(:,1),xyYFP(:,2)); %get indexes of cells within YFP pos gate
end

%% gate on mcherry cells looking at cells + single
fh = figure;
hold on
for i = 1:length(subset)
    plot(FSCA{i}(tfSingleCells{i}),mch{i}(tfSingleCells{i}),'o','MarkerSize',2,'MarkerEdgeColor',myColor(i,:),'MarkerFaceColor',myColor(i,:)) %scatter plot of YFP signal in single cells
end
yline(maxSigVal);
set(gca,'Yscale','log','Xscale','log');
xlabel('FSC-A'); ylabel('mCherry-A');
xlim([10 100000]);
grid on;
legend(description(subset),'location','eastoutside');
posmCh = drawpolygon(gca); %draw rectangle to gate mcherry cells
xymCh = posmCh.Position; %get coordinates of rectangle
close(fh)
tfmCh = {};
for i = 1:length(subset)
    tfmCh{i} = inpolygon(FSCA{i},mch{i},xymCh(:,1),xymCh(:,2)); %get indexes of cells within mcherry pos gate
end

%% sanity check: plot sum of cells in each gate (should be smaller with each subsequent gating step)
figure;
for i = 1:length(subset)
    curData = [length(FSCA{i}),sum(tfCells{i}),sum(tfSingle{i}),sum(tfSingleCells{i}),sum(tfmCh{i}&tfSingleCells{i}),sum(tfCFP{i}&tfSingleCells{i}),sum(tfYFP{i}&tfSingleCells{i})];
    subplot(3,3,i)
    bar(curData)
    grid on
    title(description(subset(i)));
    set(gca,'xticklabel',{'total events','Cells','Single','Single + Cells','mCherry + Single + Cells','CFP + Single + Cells','YFP + Single + Cells'})
    xtickangle(45)
end

%% save gating coordinates to use in the future
xyCoord.Cells = xyCells;
xyCoord.Single = xySingle;
xyCoord.CFP = xyCFP;
xyCoord.YFP = xyYFP;
xyCoord.mCherry = xymCh;
xyCoord.normYFP = xyNormYFP;

save('xyCoordinate.mat','xyCoord');
