% This script scans geonmes of E. coli strain and calculates skew in trinucleotide frequency linked to colibactin mutation signature

%% user definitions
pks_hits_cutoff = 9; % num of appearances for the "colibactin" string in a proteome to be considered a pks+ strain

%% load strain/genome/phylogroup lookup table
load dataTable_Ecoli.mat; %load datatable for all the E. coli genomes
genomeFolder = '../genomes/';

% Define all possible tri-nucleotides
tri_nucleotides = {};
nucleotides = ['A', 'C', 'G', 'T'];
for i = 1:length(nucleotides)
    for j = 1:length(nucleotides)
        for k = 1:length(nucleotides)
            tri_nucleotides{end+1} = strcat(nucleotides(i), nucleotides(j), nucleotides(k));
        end
    end
end

% get index of desired trinucleotide sequences
inx_AAA = find(strcmp('AAA',tri_nucleotides));
inx_TTT = find(strcmp('TTT',tri_nucleotides));
inx_ATA = find(strcmp('ATA',tri_nucleotides));
inx_TAT = find(strcmp('TAT',tri_nucleotides));

%% Iterate over all fastq files (contigs) and measure trinuc frequencies
categories = unique(dataTable_Ecoli.phylogroups_large); %get the phylogroups
nCategories = length(categories);

%set holders for new variables
updatedAccesions = cell(size(dataTable_Ecoli.phylogroups_large,1),1);
tfAssemblyDownloaded = false(size(dataTable_Ecoli.phylogroups_large,1),1);
tfUpdatedAssembly = false(size(dataTable_Ecoli.phylogroups_large,1),1);
ratio_ATA_AAA = nan(size(dataTable_Ecoli.phylogroups_large,1),1);

tic;
for iCat =1:nCategories
    curCat = categories{iCat};%each phylogroup at a time
    curFolder = [genomeFolder curCat '/'];%create subfolder for each phylogroup
    fastaFiles = dir([curFolder  '*.fna']); %get all files for current phylogroup

    strTitle = ['Phylogroup: ' curCat ' (' num2str(iCat) ' of ' num2str(nCategories) ')'];
   
    % Extract the sequence data
    nFiles = length(fastaFiles); %number of genomes in phylogroup
    genomeAccesionByAnalysisOrder = cell(size(fastaFiles)); % holder for genome analysis order
    n_contigs = zeros(size(fastaFiles));; % keep track of the number of contigs (holder)

    % create place holders to trinuc frequencies
    freq_AAA_TTT = nan(size(n_contigs)); 
    freq_ATA_TAT = nan(size(n_contigs));

    parfor iFile = 1:nFiles
        strMessage = ['Testing ' num2str(iFile) ' of ' num2str(nFiles) ' genomes']; %gives progress in command window
        disp(strMessage);
        curFastaFile = fastaFiles(iFile).name;
        tokens = split(curFastaFile,'_'); %extract accession
        genomeAccesion = ['GCA_' tokens{2}]; %add 'GCA' to beginning of accession number
        genomeAccesionByAnalysisOrder{iFile} = genomeAccesion; %add accession number in order of analysis
        
        fastaStruct = fastaread([curFolder curFastaFile]); %read in sequence
        frequencies = zeros(1, length(tri_nucleotides)); % holder for frequences of each trinuc sequence
        
        for iContig = 1:length(fastaStruct) % iterate over all contigs
            dna_sequence = fastaStruct(iContig).Sequence; %sequence of current contig

            % Calculate the frequency of each tri-nucleotide
            for i = 1:length(tri_nucleotides)
                frequencies(i) = frequencies(i)+count(dna_sequence, tri_nucleotides{i});
            end
        end
        n_contigs(iFile) = iContig; %store # contigs for each genome file

        % calculate the freq of aaa/ttt vs ata/tat;
        freq_AAA_TTT(iFile) = frequencies(inx_AAA)+frequencies(inx_TTT);
        freq_ATA_TAT(iFile) = frequencies(inx_ATA)+frequencies(inx_ATA);
    end % end of parfor loop

    % calculate and store the triNuc ratio and n_contigs by the order of genome in dataTable_Ecoli
    for iFile = 1:nFiles
        genomeAccesion = genomeAccesionByAnalysisOrder{iFile}; % extract accession number in order of analysis
        inxStrain = find(strncmp(genomeAccesion,dataTable_Ecoli.AssemblyAccession,13)); % looking at the first 13 chars (IGNORING version)
        
        if(~strcmp(genomeAccesion,dataTable_Ecoli.AssemblyAccession{inxStrain}))
            tfUpdatedAssembly(inxStrain) = true; % flag that a newer genome assembly file was used
        end
        tfAssemblyDownloaded(inxStrain) = true; % mark the assembly was downlaoded and analyzed
        ratio_ATA_AAA(inxStrain) = freq_ATA_TAT(iFile)./freq_AAA_TTT(iFile); %save the ratio/skew
        genomeContigNum(inxStrain) = n_contigs(iFile); % save number of contigs for each files
        updatedAccesions{inxStrain} = genomeAccesion; %save the genome accession number used for assemblies with updated accessions
    end
end
toc

inxEmpty = find(cellfun(@isempty,updatedAccesions)); %get location of failed genome downloads
updatedAccesions(inxEmpty) = {'missing'}; % flag failed geonme donwloads as missing

%% find how many proteins in each strain are annotated with "colibactin" string
n_pks_hits = nan(size(dataTable_Ecoli.AssemblyAccession)); %holder

fileID = fopen("colibactin_counts.txt", 'r');
while ~feof(fileID)  % feof tests for end of the file
    curLine = fgetl(fileID);  % Read one line from the file
    tokens = split(curLine,{'//','_',':'}); %get accession
    curAccesion = ['GCA_' tokens{3}]; 
    cur_n_pks_hits = str2num(tokens{end}); %extract number of colibactin annotated genes (end of each line in the file)
    curInx = find(ismember(updatedAccesions,curAccesion)); %see if genome was analyzed for trinuc
    n_pks_hits(curInx) = cur_n_pks_hits; % assign number of colibactin annotated proteins in order of files analyzed
end
fclose(fileID);

%% filter genomes by those that are pks+ by cutoff and plot skew with pks+ labeled
categories_reordered = categories([5 6 1 4 2 7 8 3]); % order by phylogeny (Fig 4, PMID: 33500552)
h = figure('color','white'); hold on;

% make holders
labels = {};
totalGenomeCounter = 0;
n_pks_pos = zeros(size(categories_reordered));
n_pks_neg = zeros(size(categories_reordered));
n_genomes = zeros(size(categories_reordered));

for iCat =1:length(categories_reordered)
    curCat = categories_reordered{iCat}; % select a phylogroup
    inxStrains = find(strcmp(curCat,dataTable_Ecoli.phylogroups_large)); %get index of genomes in current phylogroup
    r = ratio_ATA_AAA(inxStrains); % trinuc skew/ratio for current genomes
    cur_n_pks_hits = n_pks_hits(inxStrains); % get pks+ status for current genomes
    tf_pks = cur_n_pks_hits>pks_hits_cutoff; %logical for genomes meeting pks+ cutoff
    if(sum(tf_pks)) % some genomes are pks+
        % plot pks+
        errorbar(iCat-0.05,nanmean(r(tf_pks)),nanstd(r(tf_pks)),'or','MarkerFaceColor','r');
        plot(iCat-0.05,r(tf_pks),'.r');
        % plot pks-
        errorbar(iCat+0.05,nanmean(r(~tf_pks)),nanstd(r(~tf_pks)),'ok','MarkerFaceColor','k');
        plot(iCat+0.05,r(~tf_pks),'.k');
    else %if there are no pks+ genomes in the phylogroup
        errorbar(iCat,nanmean(r),nanstd(r),'ok','MarkerFaceColor','k');
        plot(iCat,r,'.k');
    end

    n_pks_pos(iCat) = sum(tf_pks); %save number of pks+ genomes in phylogroup
    n_pks_neg(iCat) = sum(~tf_pks); %save number of pks- genomes in phylogroup
    n_genomes(iCat) = length(tf_pks);% save total number of genomes in phylogroup
    labels{iCat} = [curCat ' (N=' num2str(sum(~isnan(r))) ')'];
    totalGenomeCounter = totalGenomeCounter+sum(~isnan(r)); %sum of all genomes across all phylogroups
end
set(gca,"XLim",[0.5 length(categories_reordered)+0.5],'xtick',[1:8],'XTickLabel',labels,'XTickLabelRotation',45);
title(['Trinuc signature by genome (N=' num2str(totalGenomeCounter) ')']);
xlabel('Phylogroup');
ylabel('Freq_A_T_A_+_T_A_T / Freq_A_A_A_+_T_T_T');
grid on; box on;
set(gcf,'position',[1000 847 848 391]);
saveas(h,'trinuc_skew_distributions.svg');

%% visualize as violin plot
%requires Violin.m and violinplot.m which can be downloaded from
%github.com/bastibe/Violinplot-Matlab

tf_valid_r = ~isnan(ratio_ATA_AAA); % extract index for ratios for analyzed genomes
cat = dataTable_Ecoli.phylogroups_large(tf_valid_r); %extract genomes with valid ratios
r = ratio_ATA_AAA(tf_valid_r)'; %extract valid ratios
r_med = median(r); %get the median ratio for all genomes

% violins for all phylogroups
h = figure; 
subplot(2,2,[1 2]); hold on; 
violinplot(r, cat,'ShowData',false,'GroupOrder',categories_reordered,'ViolinAlpha',0.5);
yline(r_med);
set(gca,'ylim',[0.77 0.85]);
grid on;

% violins for pks+ and pks- B2 phylogroup genomes
tf_B2 = strcmp(dataTable_Ecoli.phylogroups_large,'B2'); %index for B2 phylogroup genomes
r1 = ratio_ATA_AAA(tf_B2 & tf_valid_r & n_pks_hits<=pks_hits_cutoff); %pks- ratios for B2 phylogroup
r2 = ratio_ATA_AAA([tf_B2 & tf_valid_r & n_pks_hits>pks_hits_cutoff]); %pks+ ratios for B2 phylogroup
cat1 = repmat({'pks-'},size(r1)); % assign pks- for all pks- ratios
cat2 = repmat({'pks+'},size(r2)); % assign pks+ for all pks+ ratios
r = [r1;r2]; cat = {cat1{:},cat2{:}};
subplot(2,2,3); hold on; 
violinplot(r, cat,'ShowData',false,'GroupOrder',{'pks-','pks+'},'ViolinAlpha',0.5,'MarkerSize',5);
yline(r_med);
set(gca,'ylim',[0.77 0.85]);
label1 = ['B2/pks- (N=' num2str(length(r1)) ')'];
label2 = ['B2/pks+ (N=' num2str(length(r2)) ')'];
labels = {label1,label2};
grid on;

% plot histogram of B2 pks+ and pks-
subplot(2,2,4); hold on; 
histogram(r1,[0.75:0.005:0.85],'Normalization','probability');
histogram(r2,[0.75:0.005:0.85],'Normalization','probability');
legend(labels);
grid on;

set(gcf,'position',[1000 847 848 391]);
saveas(h,'trinuc_skew_violins.svg');

%% things to save in file
inxValid = find(tfAssemblyDownloaded);
infoTable = table;
infoTable.species = dataTable_Ecoli.x_Organism_Name(inxValid);
infoTable.strain = dataTable_Ecoli.Strain(inxValid);
infoTable.phylogroup = dataTable_Ecoli.phylogroups_large(inxValid);
infoTable.phylogroup_sub = dataTable_Ecoli.Phylogroup(inxValid);
infoTable.accesion = updatedAccesions(inxValid);
infoTable.genome_size = dataTable_Ecoli.Size_Mb_(inxValid);
infoTable.n_contigs = genomeContigNum(inxValid)';
infoTable.n_pks_hits = n_pks_hits(inxValid);
infoTable.is_pks_pos = num2str([infoTable.n_pks_hits>pks_hits_cutoff]);
infoTable.colibactin_skew = num2str(ratio_ATA_AAA(inxValid));
writetable(infoTable,'trinuc_skew_master_table.xlsx');

%% save useful information
save matlab_freeze 

