%% user variables
contextLength = 5; %positions upstream or downstream of mutated site
freqCutoff = 0.10; %each sample is 4 replicates pooled, so using cutoff of 0.1 for significant frequency

%% Load the mutation tables and Ref genome
% inputFile = '10Exp_WGS.xlsx'; % file with all identified mutations (generated manually)
dataStartLoc = 2;
inputFile = 'master_mutation_final.xlsx'; %read in master mutation sheet info
data = readtable(inputFile,'sheet',1,'DataRange',dataStartLoc); %all data
rawT{1} = data(contains(data.Var1,'pks+'),:); %pks+
rawT{2} = data(~contains(data.Var1,'pks+'),:); %pks-
labels = {'pks+','pks-'};

REF = getgenbank('NZ_CP009273','FileFormat','FASTA'); %get BW25113 reference genome


%% remove rows with deletions or insertions instead of point mutations and icd "mutations" since those are artifacts of the e14 deletion
for iT = 1:length(rawT)
    curT = rawT{iT};
    mutInx = [contains(curT{:,8},'icd') | ~contains(curT{:,13},'sbs')]; %find index of icd mutations and everything that is not sbs
    mutT{iT} = curT(~mutInx,:); %remove indexes from line above
end
%% remove rows with low frequency mutations
for iT = 1:length(mutT)
    curT = mutT{iT};
    freqT{iT} = curT(curT.Var5>freqCutoff,:); %take only mutations meeting the frequency cutoff
end

%% find unique mutations in each condition (because so many samples were sent per condition)
for iT = 1:length(freqT)
    curT = freqT{iT};
    vars = table(curT.Var6, curT.Var7); %pull out AA mutation and gene
    [C, ia, ic] = unique(vars,'rows'); %get unique rows for both AA mutation and gene (possible to have same AA mutation in different genes)
    T{iT} = curT(ia,:);
end

%% add columns for original and new bases and exclude deletions
old = {}; new = {};
for iT = 1:length(T)
    curT = T{iT};
    old = extractBefore(curT{:,4},'>'); % original nucleotide
    new = extractAfter(curT{:,4},'>'); % mutant nucleotide
    curT.Var14 = old;
    curT.Var15 = new;
    T{iT} = curT;
end
    
%% find uniq mutation positions by overlap in pks- and pks+ positions
% (likely existed in the ancenstor already)

pos1 = T{1}.Var3; pos2 = T{2}.Var3;  
pos_overlap = intersect(pos1,pos2); % find mutated genomic positions shared between pks+ and pks- conditions

%delete old versions of files
delete(strcat('mutationContext_',num2str(contextLength),'.txt'));
delete(strcat(labels{1},'_context_',num2str(contextLength),'.txt'));
delete(strcat(labels{2},'_context_',num2str(contextLength),'.txt'));
delete pks+_context.txt; delete pks-_context.txt; 

myfile = fopen(strcat('mutationContext_',num2str(contextLength),'.txt'),'a'); %create new file for mutation context
seqContext = {};
for iT = 1:length(T)
    curT = T{iT};
    nMutations = size(curT,1); % discard rows with NaN at the end of the table
    for iPos = 1:nMutations
        curPos = curT.Var3(iPos); % position of current mutatino
        nucOld = curT.Var14(iPos); % original nucleotide
        nucNew = curT.Var15(iPos); % mutated nucleotide
        nucOldLength = numel(nucOld); % size of mutated site
        prefixEndPos = curPos-1; % position up to the mutation site
        suffixStartPos = curPos+nucOldLength; % position right after the mutation site

        % annotate the mutation site (only if wasn't in ancestor strain)
        if(~ismember(curPos,pos_overlap)) % focus on new mutations (not in ancestor/pks-)
            prefixSeq = REF.Sequence((prefixEndPos-contextLength):(prefixEndPos)); %get genomic sequence upstream of mutation
            suffixSeq = REF.Sequence((suffixStartPos):(suffixStartPos+contextLength)); %get genomic sequence downstream of mutation
            mutationSiteInRef = REF.Sequence((prefixEndPos+1):(suffixStartPos-1)); %get reference nucleotide
            if(~strcmp(mutationSiteInRef,nucOld)), display('Error: sequence at site doesnt fit the REF genome'); end %make sure the reference nucleotide matches the breseq reference nucleotide

            
            siteOld = [lower(prefixSeq) nucOld{:} lower(suffixSeq)]; %sequence context of original/reference site
            siteNew = [lower(prefixSeq) nucNew{:} lower(suffixSeq)]; %sequence context with mutation
            
            fprintf(myfile,"%s\t%i\t%s\t%s\n",labels{iT},curPos,siteOld,siteNew); % write the variables to a file
            fastawrite(strcat(labels{iT},'_context_',num2str(contextLength),'.txt'),num2str(iPos),siteOld); %write a fasta file with the sequence context to use for motif analysis
        end
    end
end
fclose(myfile);
%% get random control sequences

ctrlPos = randi(length(REF.Sequence),10000,1); %get 10,000 random positions from 1 to the length of the genome

delete ctrl_context.txt; %delete old version
for iPos = 1:10000
    startPos = ctrlPos(iPos); %start at randomly generated position
    endPos = startPos+(2*contextLength); % get end position context length away (i.e. 2*5 = 10 bp downstream)
    seqContext = REF.Sequence(startPos:endPos); %get the sequence between these positions
    fastawrite('ctrl_context.txt',seqContext); % write this into a fasta file for motif analysis
end


