%% user definitions
%conditions = {'pks-','pks+'};

files = dir('*.czi');
conds = {};
for i = 1:length(files)
    inx(i) = contains(files(i).name,'background'); %find background file
    cur = files(i).name;
    conds{i} = string(extractBefore(cur,'-'));
end
junk = conds(inx);
n = length(junk);
conditions = replace([junk{:}],'_',' ');
files = files(inx);
%% quantify yfp/cfp in untagged colonies
bkgrnd = {};
for iFile = 1:n
    bkgrnd = [bkgrnd, quant_backgroundSignal(files(iFile).name)];
end

save bkgrnd.mat bkgrnd
%% plot yfp
% timeInx = [1:3; 4:6; 7:9]; %example index for each time
% times = {'24h','48h', '72h'}; %example times

timeInx = 1;
times = {'48h'};

yfpBgnd = []; cfpBgnd = [];
for iTime = 1:length(times) % # strains
    curTime = bkgrnd(timeInx(iTime,:));
    curConds = conditions(timeInx(iTime,:));
    % make variable holders
    allMedian = []; names = []; allCFP = []; 
    for iCond = 1:(n/length(times))
        cur = curTime{iCond};
        curName = curConds(iCond);
        allMedian = [allMedian; cur.intMedian']; %create array of yfp medians
        allCFP = [allCFP; cur.CFP']; %create array of CFP medians
        names = [names; repmat(curName,length(cur.intMedian),1)]; %create condition key for median array
        yfpBgnd(iTime,iCond) = mean(cur.intMedian);
        cfpBgnd(iTime,iCond) = mean(cur.CFP);
    end

end

%% save
% media = {'LB','M9','M9_noAA'}; %example media
median = {'M9'};
backgroundSignal = {yfpBgnd, cfpBgnd, media};
save backgroundSignal.mat backgroundSignal
