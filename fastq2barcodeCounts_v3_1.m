%% this script identifies ASKA barcode within reads in multiple fastq files.
% this script assumes both copies of the library were pooled and treated as
% one sample.
% the script can use two alternative methods:
% 1. Identify an EXACT mantch of a 20bp barcode that immedietely after a
% PS1 anchor (a conserved area upstream of the barcode). minute per screen
% 2. Serach for matches to any barcode (15-25bp) in the approx region of
% barcode. The method is slower


%% some user definitions
tfSlowMethod = true; % set to false if you want only to run with Method #1 (exact match of 20bp)

% load the barcode to gene map
load('./ASKA_lookup_map.mat'); % local copy
clear mapEven; clear mapOdd; % removing barcode lookup tables for Even and Odd barcodes (legacy)

sqCutoff = 10; % score quality cutoff
pattern_PS1 = 'AGCTGCTTCG'; % last 10 bp of PS1 (used as anchor to identify the barcode)
pattern_PS16 = 'GAATCTTCG'; % first 10 bp of PS16
%% iterate over files and extract barcode data on the fly
fastqList =dir('*_merged');

allBarcodes = map.keys; % get all barcodes
tfValidBarcode = (cellfun(@length,allBarcodes)>15 & cellfun(@length,allBarcodes)<25); % barcode length must be between 15-25bp
allValidBarcodes = allBarcodes(tfValidBarcode); %filter to barcodes of proper length


%% prealocate an empty dataset (required for parallel computing with parfor)
dataset = [];
for iFile=1:length(fastqList)
    dataset(iFile).nReads = 0;
    dataset(iFile).nPosEven = 0; dataset(iFile).nPosOdd = 0; dataset(iFile).nPosMerged = 0;
    dataset(iFile).hits = {};
    dataset(iFile).fileName = '';
    dataset(iFile).totalRuntime = 0;
    dataset(iFile).normRuntime = 0;
    dataset(iFile).datatime = datetime;
end


%% Iterate and identify barcodes in each fastq file

% a few useful variables (identical for all fastq files
uniqueGeneNames = sort(unique(ASKA.geneName));
[~,inxInASKA] = ismember(uniqueGeneNames,ASKA.geneName);

% Master loop!
parfor iFile=1:length(fastqList)
    tic;
    % Import data from the fastq file
    [~, Sequence, Qual] = fastqread(fastqList(iFile).name);
    
    % initialize empty variables to hold the data
    nReads = length(Sequence);
    
    hits = {};
    hits.geneNames = uniqueGeneNames;
    hits.mergedCounts = zeros(size(hits.geneNames)); % screen strain counts (merged for even and odd)
    hits.evenCounts = zeros(size(hits.geneNames)); % screen strain results (only even strains)
    hits.oddCounts = zeros(size(hits.geneNames)); % screen strain results (only even strains)
    
    hits.ASKA = ASKA; % hold all infromation on the ASKA library (7259 rows, by barcode)
    hits.ASKA.barcodeCounts = zeros(size(ASKA.barcode)); % add a new field to ASKA
    
    %% Iterate over all sequences in fastq file and identify barcode
    
    for i=1:nReads
        curBarcode = '';
        % get the sequence and mask low quality region
        curSeq = Sequence{1,i};
        
        if(length(curSeq)>=51) % only use reads of certain length
            curQS_str = Qual{1,i};
            curQS = double(curQS_str) - 33; % covert string to interger quality score
            tf = (curQS<sqCutoff); % logical if meets quality cutoff
            curSeqMasked = curSeq;
            curSeqMasked(tf) = 'N';
            
            % METHOD (1): identify the barcode by PS1 anchor and extract 20bp
            k = strfind(curSeq(1:30),pattern_PS1);
            curBarcode = seqrcomplement(curSeqMasked((k+10):(k+29)));
            
            if(map.isKey(curBarcode)) % barcode was easily found (by METHOD 1)
                geneName = ASKA.geneName(map(curBarcode));
                hits.ASKA.barcodeCounts(map(curBarcode)) = hits.ASKA.barcodeCounts(map(curBarcode))+1;
            elseif(tfSlowMethod) % barcode not found by METHOD 1, try METHOD 2 (if tfSlowMethod is true by user's choice)
                curBarcodeRegion = seqrcomplement(curSeqMasked(21:51)); % based on primer design the barcode should be at
                
                if(contains(curBarcodeRegion,allValidBarcodes))
                    for iBarcodeQuery=1:length(map)
                        if(tfValidBarcode(iBarcodeQuery))
                            curBarcodeQuery = allBarcodes{iBarcodeQuery};
                            k = strfind(curBarcodeRegion,curBarcodeQuery); %search for barcode in barcode region
                            if(~isempty(k)) % the curBarcodeRegion matches the curBarcodeQuery
                                queryLength = length(curBarcodeQuery);
                                curBarcode = curBarcodeQuery;
                                break;
                            end
                        end
                    end
                end
                if(map.isKey(curBarcode))
                    geneName = ASKA.geneName(map(curBarcode));
                    hits.ASKA.barcodeCounts(map(curBarcode)) = hits.ASKA.barcodeCounts(map(curBarcode))+1;
                else
                    % no match was found by Method 1 or 2
                end
            end
        else
            % no match, sequence too short
        end
        
        if(~mod(i,100000))
            fprintf('File %d of %d / seq %d of %d\n',iFile,length(fastqList),i,nReads); %status of how many reads per file have been analyzed
        end
    end
    
    %% Calculate all variables from the barcodeCounts
    
    % Counts of the barcodes
    hits.ASKA.barcodeCounts; % was just calculated
    
    % Iterate over all barcode and calculate counts per strain
    for iBarcode = 1:length(hits.ASKA.barcodeCounts)
        curBarcode = hits.ASKA.barcode(iBarcode);
        curCount = hits.ASKA.barcodeCounts(iBarcode);
        curGeneName = hits.ASKA.geneName(iBarcode);
        curPlate = hits.ASKA.plate{iBarcode};
        [~,inxGeneName] = ismember(curGeneName,hits.geneNames); % find the indexs of curGeneName [1-3680]
        if ~isnan(curPlate)
            if(mod(curPlate,2)) % if plate number is odd
                hits.oddCounts(inxGeneName) = curCount;
                hits.mergedCounts(inxGeneName) = hits.mergedCounts(inxGeneName) + curCount;
            else % plate number is even
                hits.evenCounts(inxGeneName) = curCount;
                hits.mergedCounts(inxGeneName) = hits.mergedCounts(inxGeneName) + curCount;
            end
        end
    end
    
    % Relative frequency of counts (Read Per Million)
    nPosEven = sum(hits.evenCounts); nPosOdd = sum(hits.oddCounts); nPosMerged = nPosEven+nPosOdd;
    hits.mergedRPM = hits.mergedCounts/nPosMerged*1000000;
    hits.oddRPM = hits.oddCounts/nPosOdd*1000000;
    hits.evenRPM = hits.evenCounts/nPosEven*1000000;
    
    % save the hits variable into a array of structures called dataset
    dataset(iFile).nReads = nReads;
    dataset(iFile).nPosEven = nPosEven; dataset(iFile).nPosOdd = nPosOdd; dataset(iFile).nPosMerged = nPosMerged;
    dataset(iFile).hits = hits;
    dataset(iFile).fileName = fastqList(iFile).name;
    dataset(iFile).totalRuntime = round(toc/60/60,2);
    dataset(iFile).normRuntime = dataset(iFile).totalRuntime/nReads*1000000; % run time for millon reads
    dataset(iFile).datatime = datetime;
    
    % print how long it took per million reads for this file
    fprintf('File %d of %d is done (%3.2f hr at %3.2f hr per M reads)\n',iFile,length(fastqList),dataset(iFile).totalRuntime, dataset(iFile).normRuntime);
    
end

%% save matlab structure (dataset) and export counts tables as csv files

save dataset_v3 dataset

%% export csv files

% prepare matrices that will hold the exported data
r = length(uniqueGeneNames); c = length(length(fastqList));
colLabelsEven = cell(1,c); colLabelsOdd = cell(1,c); colLabelsMerged = cell(1,c);
countsEven = zeros(r,c); countsOdd = zeros(r,c); countsMerged = zeros(r,c);

for iFile=1:length(fastqList)
    colLabelsEven{iFile} = ['Even_' fastqList(iFile).name];
    countsEven(:,iFile) = dataset(iFile).hits.evenCounts;
    colLabelsOdd{iFile} = ['Odd_' fastqList(iFile).name];
    countsOdd(:,iFile) = dataset(iFile).hits.oddCounts;
    colLabelsMerged{iFile} = ['OddPlusEven_' fastqList(iFile).name];
    countsMerged(:,iFile) = dataset(iFile).hits.mergedCounts;
end

% export the data to excel files
writematrix('Strains','countsOddPlusEven.xlsx','Range','A1');
writecell(colLabelsMerged,'countsOddPlusEven.xlsx','Range','B1');
writecell(uniqueGeneNames,'countsOddPlusEven.xlsx','Range','A2');
writematrix(countsMerged,'countsOddPlusEven.xlsx','Range','B2');

% convert excel files to csv and delete the excel files
[~,~,rawOddPlusEven] = xlsread('countsOddPlusEven.xlsx');
writecell(rawOddPlusEven, 'countsOddPlusEven.csv')
delete 'countsOddPlusEven.xlsx'


%% prepare some QA plots

% sequenicing coverage
h = figure;
RB = redbluecmap(31);
myRed = RB(2,:); myBlue = RB(10,:); myGray = [0.8 0.8 0.8];

subplot(3,1,1); hold on;
x = 1:length(fastqList);
yTotal = [dataset.nReads];
yOdd = [dataset(:).nPosOdd];
yMerged = [dataset(:).nPosMerged];
bar(x,yTotal,'FaceColor',myGray);
bar(x,yMerged,'FaceColor',myRed);
bar(x,yOdd,'FaceColor',myBlue);
legend({'Total','Even','Odd'},'location','southeast');
xlabel('Sample #'); xtickangle(90);
ylabel('# of mapped reads');
set(gca,'xtick',[1:length(yTotal)]);
grid on; box;

subplot(3,1,2); hold on;
x = 1:length(fastqList);
yTotal = ones(size(yMerged));
yOdd = [dataset(:).nPosOdd]./[dataset.nReads];
yEven = [dataset(:).nPosEven]./[dataset.nReads];
yMerged = yOdd+yEven;
bar(x,yTotal*100,'FaceColor',myGray);
bar(x,yMerged*100,'FaceColor',myRed);
bar(x,yOdd*100,'FaceColor',myBlue);
xlabel('Sample #'); xtickangle(90);
ylabel('% of mapped reads');
set(gca,'ytick',[0:10:100],'xtick',[1:length(yTotal)]);
grid on; box;

subplot(3,1,3); hold on;
x = 1:length(fastqList);
yOdd = [dataset(:).nPosOdd];
yEven = [dataset(:).nPosEven];
yMerged = yOdd+yEven;
yLFR = log2(yOdd./yEven);
stem(x,yLFR,'color','k','markerfacecolor','yellow');
xlabel('Sample #'); xtickangle(90);
ylabel('odd/even ratio (log_2)');
myLim = 1.2*max(abs(yLFR))
set(gca,'ylim',[-myLim myLim],'xlim',[0 max(x)+1]),'xtick',[1:length(yLFR)];
set(gca,'xtick',[1:length(yTotal)]);
grid on; box;

set(gcf,'position',[20         150        1550         593]);
saveas(h, 'Coverage statistics', 'png');

% run/computation efficeny
h = figure
subplot(2,1,1); hold on;
y = [dataset.totalRuntime];
bar(y,'facecolor',myGray);
xlabel('Sample #'); xtickangle(90);
ylabel('Total run time (hr)');
set(gca,'xtick',[1:length(y)]);
grid on; box on;

subplot(2,1,2); hold on;
y = [dataset.normRuntime];
bar(y,'facecolor',myGray);
xlabel('Sample #'); xtickangle(90);
ylabel('Barcode mapping rate (hr per 10^6 reads)');
set(gca,'xtick',[1:length(y)]);
set(gcf,'position',[70         100        1550         593]);
grid on; box on;

saveas(h, 'Mapping efficency', 'png');


