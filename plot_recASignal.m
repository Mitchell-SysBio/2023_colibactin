%% user defined variables
nTime = 1; %number of time points
times = {'24h'};
noTouchCutOff = 10; % distance in microns to say colonies are not touching 
touchCutoff = 4; %max distance in microns between colnoies to say they are touching
condsOrder = {'EcN','BAC'};
%% read in files
files = dir('*.czi');

for i = 1:length(files)
    fileName{i} = files(i).name;
end

%% segment Nissle strain microscopy images
EcNPos = {}; EcNNeg = {};
for i = 1:length(files)
    if contains(fileName{i},'EcN') & contains(fileName{i},'pos')
        EcNPos = [EcNPos, segment_recASignal(fileName{i})]; % apply recA function to pks+ samples
    elseif contains(fileName{i},'EcN') & contains(fileName{i},'neg')
        EcNNeg = [EcNNeg, segment_recASignal(fileName{i})]; %apply recA function to pks- samples
    end
end

EcN = {EcNPos; EcNNeg};

%% segment BAC strain microscopy images
BACPos = {}; BACNeg = {};
for i = 1:length(files)
    if contains(fileName{i},'BAC') & contains(fileName{i},'pos')
        BACPos = [BACPos, segment_recASignal(fileName{i})]; %apply recA function to pks+ samples
    elseif contains(fileName{i},'BAC') & contains(fileName{i},'neg')
        BACNeg = [BACNeg, segment_recASignal(fileName{i})]; %apply recA function to pks- samples
    end
end

BAC = {BACPos; BACNeg};
%% save conds
conds = {EcN; BAC};
save EcN.mat EcN; %save the EcN conds variable 
save BAC.mat BAC; %save the BAC conds variable
save conds.mat conds; %save the entire conds variable

%% extract trace from peak YFP signal of touching and non touching for both pos and neg
% get all traces and means with error in subplots for non touching pks+,
% touching pks+, touching pks-
decayMat = []; mat = []; nonTouch = [];
allPosPos = {}; allYFPPos = {}; allRepPos = {};
allPosNeg = {}; allYFPNeg = {}; allRepNeg = {};
allDecay = []; allDist = [];

for iStrain = 1:2 %loop through each strain
    curStrain = conds{iStrain};
    for iCond = 1:2 %loop through pks+ and then pks-
        col = curStrain{iCond};
        decay = [];
        for i=1:length(col)
            scaleFactor = double(col{i}.micron_per_pixel); % microns per pixel from function
            % calculate the distance jumps between intensity points
            if ~isempty(col{i}.x)
                dx = col{i}.x(2)-col{i}.x(1); dy = col{i}.y(2)-col{i}.y(1);
                dist_unit = sqrt(dx^2+dy^2); %number of pixels between each intensity point
                % calculate location of colony edges along line
                if size(col{i}.edges,1) == 2 %exclude anything messed up by human error (i.e. not drawing a line)
                    dxRep = col{i}.x(1) - col{i}.edges(1,1); %dx of line start in reporter colony to reporter colony edge
                    dyRep = col{i}.y(1) - col{i}.edges(1,2); %dy of line start in reporter colony to reporter colony edge
                    dxTox = col{i}.x(1) - col{i}.edges(2,1); %dx of line start in reporter colony to toxic colony edge
                    dyTox = col{i}.y(1) - col{i}.edges(2,2); %dy of line start in reporter colony to toxic colony edge
                    distRep = sqrt(dxRep^2 + dyRep^2); % number of pixels between line start in reporter colony and edge of reporter colony
                    distTox = sqrt(dxTox^2 + dyTox^2); % number of pixels between line start in reporter colony and edge of toxic colony
                end
                inx = round(distTox,0);
                decay{i}.pos = [inx:-1:1]*dist_unit*scaleFactor; %get positions in um from edge of toxic colony to end of roi line
                decay{i}.rep = decay{i}.pos(1) - distRep*scaleFactor; %get distance between colonies
                decay{i}.yfp = col{i}.c3(1:inx); %extract the yfp intensity for the reporter colony
                decay{i}.edge = inx; %save index of edge of toxic colony
                decay{i}.repEdge = round(distRep,0); %save index of edge of reporter colony
            else %removes tiles where i made mistakes with drawing the rois or where the field of view was bad
                decay{i}.pos = nan;
                decay{i}.rep = nan;
                decay{i}.yfp = nan;
                decay{i}.edge = nan;
            end
        end
        
        % separate positive and negative along with contact and no contact
        pos = {}; yfp = {};
        noSignalPos = {}; noSignalYFP = {}; noSignalRep = {};
        noTouchTraceMax = {}; noTouchTrace = {}; touchTrace = {};
        noTouchDist = []; touchDist = []; 
        if iCond == 1 % pks+
            for i=1:length(col)
                % non-touching colonies
                if decay{i}.rep >= noTouchCutOff %colonies separated by more than noTouchCutOff
                    pos = [pos;decay{i}.pos]; % create cell with all position arrays
                    yfp = [yfp;decay{i}.yfp]; % create cell with all yfp arrays
                    maxYFP = find(decay{i}.yfp == max(decay{i}.yfp)); % find the index of the peak yfp signal along the profile
                    junk = decay{i}.yfp(1:decay{i}.edge); % get yfp decay from edge of tox colony
                    junk(maxYFP:decay{i}.edge) = nan; %replace values between edge of tox and peak YFP with nan
                    noTouchTrace{end+1} = flip(junk); % flip so the values go from peak yfp to end of roi line in the reporter colony
                    noTouchDist(end+1) = decay{i}.rep; %save distance between colonies in an array
                % touching colonies
                elseif decay{i}.rep <= touchCutoff %colonies separated by less than the touchCutoff
                    pos = [pos;decay{i}.pos]; % cell with all position arrays
                    yfp = [yfp;decay{i}.yfp]; %cell with all yfp arrays
                    maxYFP = find(decay{i}.yfp == max(decay{i}.yfp)); %find the index of the peak yfp signal
                    junk = decay{i}.yfp(1:decay{i}.edge); %yfp decay from edge of tox colony
                    junk(maxYFP:decay{i}.edge) = nan; %replace values from the edge to the max YFP with nan
                    touchTrace{end+1} = flip(junk); % flip so the values go from peak yfp to end of roi line in the reporter colony
                    touchDist(end+1) = decay{i}.rep; %save distance between colonies in an array
                end
                allPosPos = [allPosPos; decay{i}.pos]; %create cell for all pks+ position arrays
                allYFPPos = [allYFPPos; decay{i}.yfp]; %create cell for all pks+ yfp arrays
                allRepPos = [allRepPos; decay{i}.rep]; %create cell for all pks+ distances between colonies
            end
        else % pks- (cannot take from max yfp signal because the signal is too flat across entire colony)
            for i=1:length(col)
                %touching colonies only
                if decay{i}.rep <= touchCutoff % colonies separated by less than the touchCutoff
                    pos = [pos;decay{i}.pos]; % cell with all position arrays
                    yfp = [yfp;decay{i}.yfp]; %cell with all yfp arrays
                    junk = decay{i}.yfp(1:decay{i}.edge); %yfp decay from edge of tox colony
                    junk(decay{i}.repEdge:decay{i}.edge) = nan; %replace values from the edge of tox to edge of reporter with nan
                    TouchTrace{end+1} = flip(junk); %flip to get yfp values from edge of reporter colony to end of line inside reporter colony
                    TouchDist(end+1) = decay{i}.rep; %save the distance between colonies in an array
                %still save non touching colonies in case
                elseif decay{i}.rep >= noTouchCutOff %colonies separated by more than noTouchCutOff
                    pos = [pos;decay{i}.pos]; % cell with all position arrays
                    yfp = [yfp;decay{i}.yfp]; %cell with all yfp arrays
                    junk = decay{i}.yfp(1:decay{i}.edge); %yfp decay from edge of tox colony
                    junk(decay{i}.repEdge:decay{i}.edge) = nan; %replace values from the edge of tox to edge of reporter with nan
                    noTouchTrace{end+1} = flip(junk); %flip to get yfp values from edge of reporter colony to end of line inside reporter colony
                    nTouchDist(end+1) = decay{i}.rep; %save the distance between colonies in an array
                end
                allPosNeg = [allPosNeg; decay{i}.pos]; %create cell for all pks- position arrays
                allYFPNeg = [allYFPNeg; decay{i}.yfp]; %create cell for all pks- yfp arrays
                allRepNeg = [allRepNeg; decay{i}.rep]; %create cell for all pks- distances between colonies
            end
        end
        decayMat{iCond} = {noTouchTrace; touchTrace}; %combine non-touching and touching yfp traces to one variable
        distMat{iCond} = {noTouchDist; touchDist}; %combine non-touching and touching colony distances to one variable
    end
    allDecay{iStrain} = decayMat; % save yfp decays for each strain
    allDist{iStrain} = distMat; % save colony distances for each strain
end
%% plot
titles = {'non-touching pks+','touching pks+','non-touching pks-','touching pks-'};
distCutoff = 150; %minimum length of roi lines (measured points, not microns)
ymax = [0.1, 0.2]; %upper y lim by strain

for iStrain = 1:2
    decayMat = allDecay{iStrain};
    figure; hold on;
    for iCond = 1:2
        %set up subplot indexing by condition
        for iTouch = 1:2 %1 is no touch, 2 is touch
            if iCond ==1 & iTouch == 1 %no touch pks+
                inx = 1;
            elseif iCond == 1 & iTouch == 2 %touch pks+
                inx = 2;
            elseif iCond ==2 & iTouch ==1 %no touch pks-
                inx = 3;
            else
                inx = 4; %touch pks-
            end
            curMat = []; normMat = [];
            subplot(2,4,inx)
            for i = 1:length(decayMat{iCond}{iTouch})
                cur = decayMat{iCond}{iTouch}{i}; % extract yfp decay for each cell by each condition
                if length(cur) < distCutoff %find lines shorter than the distance we want to plot
                    cur(end+1:distCutoff) = nan; %add nans to end of data to make it the distance we want to plot
                else
                    cur = cur(1:distCutoff); %make sure all are the same length by removing values after distance cutoff
                end
                curMat = [curMat; cur]; %put all rois with nans into matrix
            end
            meanMat = nanmean(curMat,1); %get mean (excluding nans)
            stdMat = nanstd(curMat,[],1); %get standard deviation
            halfDecay = (max(meanMat)-nanmean(meanMat((distCutoff-20):distCutoff)))/2; %find what half the max signal is
            halfInx = max(find(meanMat >= (halfDecay+nanmean(meanMat((distCutoff-20):distCutoff))))); %index where half max signal is found (distance)
            x = 1:scaleFactor:(distCutoff*scaleFactor); %set x variable in microns
            % plot traces of all replicates and mean overlaid 
            plot(x,curMat,'b'); hold on; %plot all traces
            plot(x,meanMat,'k','LineWidth',3) %plot mean of traces
            grid on; box on;
            title(titles{inx})
            xlim([1 distCutoff*scaleFactor]);
            ylim([0 ymax(iStrain)]);
            text(100,ymax(iStrain)-0.02,['n = ' num2str(size(curMat,1))]) % label number of colonies analyzed
    
            subplot(2,4,(inx+4)) %plot shaded error bars
            shadedErrorBar(x,meanMat,stdMat)
            grid on; box on;
            title(titles{inx})
            xlim([1 distCutoff*scaleFactor]);
            ylim([0 ymax(iStrain)]);
            xline(halfInx*scaleFactor); % mark 50% of the maximum effect
        end
    end
    sgtitle(condsOrder{iStrain})
end

%% t-test touching vs non-touching at different distances for EcN
x100 = []; x150 = []; x200 = []; x300 = [];
y100 = []; y150 = []; y200 = []; y300 = [];

for iStrain = 1:2
    decayMat = allDecay{iStrain};
    for i = 1:length(decayMat{1}{1}) %non-touching colonies
        cur = decayMat{1}{1}{i};
        x100(i) = cur(25); %100um index
        x150(i) = cur(41); %150um index
        x200(i) = cur(55); %200um index
        x300(i) = cur(83); %300um index
    end
    
    for i = 1:length(decayMat{1}{2}) %touching colonies
        cur = decayMat{1}{2}{i};
        y100(i) = cur(25); %100um index
        y150(i) = cur(41); %150um index
        y200(i) = cur(55); %200um index
        y300(i) = cur(83); %300um index
    end
    
    figure; 
    subplot(4,1,1); hold on; 
    % plot overlaid histograms of 100um yfp intensity and # of colonies with signal at that distance
    histogram(x100(~isnan(x100)),[0:0.005:0.07],'Normalization','probability')
    histogram(y100(~isnan(y100)),[0:0.005:0.07],'Normalization','probability')
    text(0.06,0.1,num2str(sum(~isnan(x100))),'Color','b')
    text(0.06,0.2,num2str(sum(~isnan(y100))),'Color','r')
    box on; grid on;
    title([condsOrder(iStrain), ' 100um'])
    
    subplot(4,1,2); hold on;
    % plot overlaid histograms of 150um yfp intensity and # of colonies with signal at that distance
    histogram(x150(~isnan(x150)),[0:0.005:0.07],'Normalization','probability')
    histogram(y150(~isnan(y150)),[0:0.005:0.07],'Normalization','probability')
    text(0.06,0.1,num2str(sum(~isnan(x150))),'Color','b')
    text(0.06,0.2,num2str(sum(~isnan(y150))),'Color','r')
    box on; grid on;
    title([condsOrder(iStrain), ' 150um'])
    
    subplot(4,1,3); hold on;
    % plot overlaid histograms of 200um yfp intensity and # of colonies with signal at that distance
    histogram(x200(~isnan(x200)),[0:0.005:0.07],'Normalization','probability')
    histogram(y200(~isnan(y200)),[0:0.005:0.07],'Normalization','probability')
    text(0.06,0.1,num2str(sum(~isnan(x200))),'Color','b')
    text(0.06,0.2,num2str(sum(~isnan(y200))),'Color','r')
    box on; grid on;
    title([condsOrder(iStrain), ' 200um'])
    
    subplot(4,1,4); hold on;
    % plot overlaid histograms of 300um yfp intensity and # of colonies with signal at that distance
    histogram(x300(~isnan(x300)),[0:0.005:0.07],'Normalization','probability')
    histogram(y300(~isnan(y300)),[0:0.005:0.07],'Normalization','probability')
    text(0.06,0.1,num2str(sum(~isnan(x300))),'Color','b')
    text(0.06,0.2,num2str(sum(~isnan(y300))),'Color','r')
    box on; grid on;
    title([condsOrder(iStrain), ' 300um'])
    legend({'non-touching','touching'})
    xlabel('YFP intensity');
    ylabel('proprotion')
    
    %t-tests at each distance per strain
    [h100(iStrain) p100(iStrain)] = ttest2(x100,y100);
    [h150(iStrain) p150(iStrain)] = ttest2(x150,y150);
    [h200(iStrain) p200(iStrain)] = ttest2(x200,y200);
    [h300(iStrain) p300(iStrain)] = ttest2(x300,y300);
end