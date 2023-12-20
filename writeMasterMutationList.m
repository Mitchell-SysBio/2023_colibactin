data = readcell('breseq_mutation_list.xlsx'); %manually created excel file of breseq mutations
REF = getgenbank('NZ_CP009273','FileFormat','FASTA'); %reference sequence for sequence near mutation
contextLength = 5; %bases to go upstream and downstream of mutation for context
%% add header names
data{1,8} = 'gene';
data{1,9} = 'synonomous';
data{1,10} = 'upstream_sequence';
data{1,11} = 'downstream_sequence';
data{1,12} = 'sbs_mutation_context';
data{1,13} = 'mutation_type';
%% clean up gene column
for i = 2:length(data)
    curMut = data{i,6};
    curGene = data{i,7};
    if ismissing(curMut) % missing values for deletions, replace with nan
        curMut = 'nan';
    end
    if contains(curMut, 'intergenic') %add gene name to new column or nan if it's an intergenic region
        data{i,8} = nan;
    else data{i,8} = curGene;
    end
end

%% add column identifying if mutations are synonymous
for i = 2:length(data)
    curData = data{i,6};
    if ismissing(curData)
        curMut = 'nan';
        data{i,6} = 'nan'; %replace missing values with nan
    else curMut = extractBefore(curData,'('); %extract ther amino acid change
    end
    if contains(curMut, 'intergenic')
        data{i,9} = nan; %label intergenic mutatations as nan
    elseif contains(curMut,'nan') %pre-labeled mutations stay nan
        data{i,9} = nan;
    elseif strcmp(curMut(1),curMut(end-1))
        data{i,9} = 'synonymous'; %if amino acid is the same (first and last letters) synonymous
    else data{i,9} = 'nonsynonymous'; %else it's nonsynonymous
    end
end

%% add sequence before and after mutation and mutation type
error = [];
for i = 2:length(data)
    curMut = data(i,:);
    curPos = curMut{3}; %mutation position
    if contains(curMut{4},'Δ')
        junk = extractAfter(curMut{4},'Δ'); %extract string after delta sign (nbp)
        junk2 = extractBefore(junk,' bp'); %extract number before bp)
        nucOldLength = str2num(junk2);  %number of nt in deletion
        if nucOldLength < 50
            data{i,13} = 'small_deletion'; %assign label of small deletion
        elseif nucOldLength > 50 & nucOldLength < 150
            data{i,13} = 'medium_deletion'; %assign label of medium deletion
        else
            data{i,13} = 'large_deletion'; %assign label of large deletion
        end
    elseif contains(curMut{4},'+') %insertions
        junk = extractAfter(curMut{4},' +');
        junk2 = extractBefore(junk,' bp');
        nucOldLength = str2num(junk2); %get number of nt in insertion
        if nucOldLength < 50
            data{i,13} = 'small_insertion'; %assign small insertion label
        elseif nucOldLength > 50 & nucOldLength < 150
            data{i,13} = 'medium_insertion'; %assign medium insertion label
        else
            data{i,13} = 'large_insertion'; %assign large insertion label
        end
    else nucOldLength = 1;
        nucOld = extractBefore(curMut{4},'>'); %extract parental bp
        nucNew = extractAfter(curMut{4},'>'); %extract mutant bp
        data{i,13} = 'sbs'; %label as sbs
    end
    prefixEndPos = curPos-1; % position up to the mutation site
    suffixStartPos = curPos+nucOldLength; % position right after the mutation site
    prefixSeq = REF.Sequence((prefixEndPos-contextLength):(prefixEndPos)); %upstream sequence
    suffixSeq = REF.Sequence((suffixStartPos):(suffixStartPos+contextLength)); %downstream sequence
    mutationSiteInRef = REF.Sequence((prefixEndPos+1):(suffixStartPos-1)); %mutated site in reference sequence
    if(~strcmp(mutationSiteInRef,nucOld))
        display('Error: sequence at site doesnt fit the REF genome');
        error = [error, i];
    end
    if nucOldLength > 2 %nan for anything that is not a single base substitution
        data{i,10} = nan;
        data{i,11} = nan;
    else
        data{i,10} = prefixSeq; %column of bp upstream of mutant base
        data{i,11} = suffixSeq; %column of bp downstream of mutant base
    end
    if nucOldLength == 1
        data{i,12} = [lower(prefixSeq) nucOld lower(suffixSeq)]; % mutation sequence context with upstream, downstream and mutated base (parental)
    else data{i,12} = nan;
    end
end


%% save new master file
writecell(data,'master_mutation_final.xlsx');

