%% user definitions
load backgroundSignal.mat; % autofluorescence of YFP and CFP at 24, 48 and 72 hours

files = dir('*.czi');
conds = {};
for i = 1:length(files)
    cur = files(i).name;
    conds{i} = string(extractBefore(cur,'-')); % assign file names to conditions
end
n = length(conds);
conditions = replace([conds{:}],'_',' '); %remove underscores from names since matlab doesn't like them


%% quantify yfp in colonies
yfp = {};
for iFile = 1:n
    yfp = [yfp, quant_recASignal(files(iFile).name)]; % (takes a few minutes per .czi file to run)
end

save yfp.mat yfp % save yfp for future use
%% plot yfp
junk = [3,1,5,7,9]; % set sample order of first time point
timeInx = [junk; junk+1]; % get indexes of samples by time point
times = {'24h','48h'};

figure; hold on;
tiledlayout(length(times),3)
nColonies = []; subYFPMat = {}; subCFPMat = {}; normYFPMat = {};
for iTime = 1:length(times) % # strains
    curTime = yfp(timeInx(iTime,:));
    curConds = conditions(timeInx(iTime,:));
    curYfpBG = backgroundSignal{1}(iTime,:); % background YFP at that time
    curCfpBG = backgroundSignal{2}(iTime,:); % background CFP at that time
    % make variable holders
    allMedian = []; names = []; allCFP = [];
    normYFP = []; subYFP = []; subCFP = []; 
    for iCond = 1:(n/length(times))
        cur = curTime{iCond};
        curName = curConds(iCond);
        if contains(curName,'LB') %set up indexes to subtract background fluorescence
            BG = 1;
        elseif contains(curName,'noAA')
            BG = 3;
        else
            BG = 2;
        end
        allMedian = [allMedian; cur.intMean']; %create array of yfp medians
        allCFP = [allCFP; cur.CFP']; %create array of CFP medians
        names = [names; repmat(curName,length(cur.intMean),1)]; %create condition key for median array
        subYFP = [subYFP; (cur.intMean-curYfpBG(BG))']; %subtract yfp background
        subYFPMat{iTime,iCond} = cur.intMean-curYfpBG(BG);
        subCFP = [subCFP; (cur.CFP-curCfpBG(BG))']; %subtract CFP background
        subCFPMat{iTime,iCond} = cur.CFP-curCfpBG(BG);
        normYFPMat{iTime,iCond} = (cur.intMean-curYfpBG(BG))./(cur.CFP-curCfpBG(BG));
        nColonies(iTime,iCond) = length(cur.CFP); % save number of colonies per condition
    end
    normYFP = subYFP./subCFP; %normalize YFP to CFP
    
    %plot YFP
    nexttile; hold on; 
    violinplot(allMedian,names,'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds));
    box on; grid on;
    title(['YFP ' times{iTime}])
    
    %plot CFP
    nexttile; hold on;
    violinplot(allCFP,names,'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds))
    box on; grid on;
    title(['CFP ' times{iTime}])
    
    %plot normalized YFP/CFP
    nexttile; hold on;
    violinplot(normYFP,names,'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds))
    box on; grid on;
    title(['normalized YFP ' times{iTime}])
    ylim([0 0.55])

end

%% make plot for paper figure
BACinx = contains(names,'BAC');
EcNInx = contains(names,'EcN');
clbNinx = contains(names,'clbN');
clbSinx = contains(names,'clbS');

% plot BAC-pks and BAC-empty
figure; hold on;
subplot(1,3,1)
violinplot(normYFP(BACinx),names(BACinx),'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds))
grid on; box on;
title('pks - ctrl')
ylim([0 0.5])

% plot EcN and clbN
subplot(1,3,3)
violinplot(normYFP(EcNInx & ~clbSinx),names(EcNInx & ~clbSinx),'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds(3:4)))
grid on; box on;
title('WT - clbN')
ylim([0 0.5])

% plot EcN and clbS
subplot(1,3,2)
violinplot(normYFP(EcNInx & ~clbNinx),names(EcNInx & ~clbNinx),'ViolinAlpha',0.3,'ShowBox',false,'ShowWhiskers',false,'GroupOrder',cellstr(curConds([3,5])))
grid on; box on;
title('WT - clbS')
ylim([0 0.5])

%% conduct t-tests
colInx = [2,1; 3,4; 3,5]; % indexes of conditions to run t-test on

h = []; p = [];
for iComp = 1:length(colInx)
    curComp = normYFPMat(2,colInx(iComp,:));
    [h(iComp) p(iComp)] = ttest2(curComp{1},curComp{2});
end

