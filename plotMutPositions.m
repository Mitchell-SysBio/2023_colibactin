%% Definitions
genome_length = 4631445;

features = readcell('genomeFeatures.xlsx'); %positions of prophages and ori/ter
macrodomains = readcell('macrodomain_positions.xlsx'); %macrodomain loci
motifInx = readmatrix("motif1_indexes.txt"); %indexes of mutations in A/T motif

T = readcell('mutationContext_5.txt');

all.strains = {T{:,1}};
all.positions = [T{:,2}];
all.seqFrom = {T{:,3}};
all.seqTo = {T{:,4}};
all.inxPks = find(strcmp(all.strains,'pks+'));

positions = all.positions(all.inxPks);
positionsMotif = positions(motifInx(1:67));


%% Calculate the corresponding angles for each position
genome_length = 4631445;
angles = (positions/genome_length)*2*pi;

% Create a unit circle
theta = linspace(0, 2*pi, 1000);
y = cos(theta);
x = sin(theta);

% Plot the unit circle
figure; hold on;
plot(x, y);
axis equal;

% Mark the positions of all mutations on the circle
for i = 1:length(positions)
    plot(sin(angles(i)), cos(angles(i)), '.k', 'MarkerSize', 10);
end

% Plot the unit circle
figure; hold on;
plot(x, y);
axis equal;

% mark positions of mutations not in motif on the circle
cur = (positions(~ismember(positions,positionsMotif))/genome_length)*2*pi;
for i = 1:length(positionsMotif)
    plot(sin(cur(i)), cos(cur(i)), 'o','MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerSize',10);
end

% Plot the unit circle
figure; hold on;
plot(x, y);
axis equal;

% mark the positions of mutations occurring in the motif on the circle
myColor = [0.8 0.3 1];
cur = (positionsMotif/genome_length)*2*pi;
for i = 1:length(positionsMotif)
    plot(sin(cur(i)), cos(cur(i)), 'o','MarkerFaceColor',myColor,'MarkerEdgeColor',myColor, 'MarkerSize',10);
end


% mark position of interest (ori and ter)
% for i = 1:length(features)-1
%     keyLabels{i} = features{i+1,1};
%     keyAngles(i) = (features{i+1,4}/genome_length)*2*pi;
% end

% keyLabels = {'ORI','terA','terB','terC','terD'};
% keyAngles(1) = (pos.ori/genome_length)*2*pi;
% keyAngles(2) = (pos.terA/genome_length)*2*pi;
% keyAngles(3) = (pos.terB/genome_length)*2*pi;
% keyAngles(4) = (pos.terC/genome_length)*2*pi;
% keyAngles(5) = (pos.terD/genome_length)*2*pi;


% for i = 1:length(keyLabels)
%     plot(sin(keyAngles(i)), cos(keyAngles(i)), '.r', 'MarkerSize', 25);
%     text(sin(keyAngles(i))+0.03, cos(keyAngles(i)), keyLabels{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 10);
% end


%% frequency of mutations in each macrodomain
domainFreq = [];
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domainFreq(i) = sum(positions >= curDomain{2} & positions <= curDomain{3});
end
domainFreq(7) = domainFreq(1) + domainFreq(7);
domainFreq(1) = [];
figure;
bar(domainFreq)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations');
title('# of mutations per macrodomain')

% normalize # mutations to length of macrodomain
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domLength(i) = curDomain{3} - curDomain{2};
end

domLength(6) = domLength(1) + domLength(6); %origin domain split into 2 b/c spans beginning and end of chromosome
domLength(1) = [];

figure;
bar(domainFreq./domLength)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations/length of domain');
title('# of mutations per macrodomain')
ylim([0 7e-5]);

%% frequency of mutations in enriched motif in each macrodomain
domainFreq = [];
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domainFreq(i) = sum(positionsMotif >= curDomain{2} & positionsMotif <= curDomain{3}); %using motif positions
end
domainFreq(7) = domainFreq(1) + domainFreq(7);
domainFreq(1) = [];
figure;
bar(domainFreq)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations');
title('# of mutations per macrodomain')

% normalize # mutations to length of macrodomain
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domLength(i) = curDomain{3} - curDomain{2};
end

domLength(6) = domLength(1) + domLength(6);
domLength(1) = [];

figure;
bar(domainFreq./domLength)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations/length of domain');
title('# of mutations per macrodomain')
ylim([0 3e-5])

%% frequency of mutations in each macrodomain excluding mutations in motif
noMotifPos = positions(~ismember(positions,positionsMotif)); % positions not in motif

domainFreq = [];
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domainFreq(i) = sum(noMotifPos >= curDomain{2} & noMotifPos <= curDomain{3});
end
domainFreq(7) = domainFreq(1) + domainFreq(7);
domainFreq(1) = [];
figure;
bar(domainFreq)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations');
title('# of non-motif mutations per macrodomain')

% normalize # mutations to length of macrodomain
for i = 1:length(macrodomains)-1
    curDomain = macrodomains(i+1,:);
    domLength(i) = curDomain{3} - curDomain{2};
end

domLength(6) = domLength(1) + domLength(6);
domLength(1) = [];

figure;
bar(domainFreq./domLength)
set(gca,'xticklabel',macrodomains(3:end,1))
box on; grid on;
xlabel('macrodomain');
ylabel('# of mutations/length of domain');
title('# of non-motif mutations per macrodomain')
ylim([0 3.5e-5]);
