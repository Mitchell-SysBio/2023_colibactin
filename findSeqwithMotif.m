%% find sequences with the enriched motif and the position of the mutated base
allData = readcell('fimo.xlsx'); %fimo results
data = allData(2:end,:);

pksContext = readcell('pks+_context_5.txt'); % fasta format of pks+ mutation context
seqName = cell2mat(pksContext(1:2:320,2));
seqContext = pksContext(2:2:320,1);
%% extract streme1 samples (most enriched motif)
seq = [];
for i = 1:189
    cur = data{i,2};
    if contains(cur,'STREME-1')
        seq = [seq,data{i,3}];
    end
end

seq = unique(seq); %some sequences fit motif multiple ways

%% save file
writecell(motifPos,'motif_sequence_context.xlsx')
writematrix(seq,'motif1_indexes.txt')