%% read in mutation data
data = readcell('master_mutation_final.xlsx');
data = data(2:end,:);
%% separate mutations by condition
pksPos = data(contains(data(:,1),'pks+'),:);
pksNeg = data(contains(data(:,1),'pks-'),:);

%% set current condition manually by user
cur = pksPos;
% cur = pksNeg;

%% extract relevant data and determine frequency of mutation types
sigData = cur(([cur{:,5}] >= 0.1),:); % meets threshold of proportion
noIcd = sigData(~contains(sigData(:,7),'icd '),:); % pull out non-icd annotated mutations
icd = sigData(contains(sigData(:,7),'icd '),:); % pull out icd annotated mutations (are actually e14 prophate)

sbs = []; small = []; medium = []; large = [];
sbs = sum(contains(noIcd(:,13),'sbs')); % get sbs number
small = sum(contains(noIcd(:,13),'small')); % small indel
medium = sum(contains(noIcd(:,13),'medium')); % medium indel
large = sum(contains(noIcd(:,13),'large')); %large indel

%% get unique deletions
uniq = []; inx = []; mutCon = {};
for i = 1:length(noIcd)
    mutCon{i} = string(noIcd{i,7});
end
[uniq inx] = unique([mutCon{:}],'stable');

%% get number of occurrence for the unique mutations
uniqMut = noIcd(inx,:);

uniqsbs = []; uniqsmall = []; uniqmedium = []; uniqlarge = [];
uniqsbs = sum(contains(uniqMut(:,13),'sbs'));
uniqsmall = sum(contains(uniqMut(:,13),'small'));
uniqmedium = sum(contains(uniqMut(:,13),'medium'));
uniqlarge = sum(contains(uniqMut(:,13),'large'));

%% get total instances of deletions as number of replicates 
x = [];
for iMut = 1:length(noIcd)
    curMut = noIcd(iMut,:);
    curFreq = curMut{5};
    if curFreq <= 0.375
        x(iMut) = 1;
    elseif curFreq > 0.375 & curFreq <= 0.625
        x(iMut) = 2;
    elseif curFreq > 0.625 & curFreq <= 0.875
        x(iMut) = 3;
    else
        x(iMut) = 4;
    end
end

if sum(x) > length(noIcd)
    print('all mutations occur in a single replicate')
else
    print('some mutations are in multiple replicates - extract those indexes for downstream analysis')
end
%% get number of replicates with icd deletion
[uniq inx] = unique(icd(:,2),'stable');
uniqIcd = icd(inx,:);
x = [];
for iMut = 1:size(uniqIcd,1)
    curMut = uniqIcd(iMut,:);
    curFreq = curMut{5};
    if curFreq <= 0.375
        x(iMut) = 1;
    elseif curFreq > 0.375 & curFreq <= 0.625
        x(iMut) = 2;
    elseif curFreq > 0.625 & curFreq <= 0.875
        x(iMut) = 3;
    else
        x(iMut) = 4;
    end
end
totalSamp = sum(x);
