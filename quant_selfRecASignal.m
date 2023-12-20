function [fluor] = quant_recASignal(fileName,tfDebug)

if(nargin == 2)
    iSampleImage = 1; % image to plot for debug purposes
else
    tfDebug = false;
end

areaMin = 10000; % to get colonies
areaMax = 175000; % to exclude masks with multiple colonies
circCutoff = 0.4; % threshold for circularity to further select only single colonies

%% open file
data = {}; 
data = bfopen(fileName); % Import the czi file (utility function from bfmatlab)
n = size(data,1);
C1 = {}; C2 = {}; C3 = {}; bwC3 = {};
for i = 1:n
    C1{i} = double(data{i,1}{1})./2^16; % brightfield
    C2{i} = double(data{i,1}{2})./2^16; % YFP
    C3{i} = double(data{i,1}{3})./2^16; % CFP
    bwC3{i} = mat2gray(data{i,1}{3}); %min/max normalize CFP and convert to double
end

%% iterate over images and make a BW mask by CFP
BW = {}; PERIM = {};
for i=1:n
    img = C3{i};
    level = graythresh(C3{i}); % get a single threshold automatically
    
    % make a bw mask by the threshold
    bw = imbinarize(img,level);
    noBorder = imclearborder(bw); %remove masked objects touching image edges
    perim = bwperim(noBorder); % bwperim find the perimeter of the obects in the bw mask
    BW{i} = noBorder;
    PERIM{i} = perim;
    
    % get properties of the bw image (identifies objects in the bw mask)
    prop = regionprops(noBorder,'Area','Centroid','PixelIdxList','Circularity');
    Props{i} = prop;
end

if tfDebug % show masks compared to original image brightfield
    figure; hold on;
    subplot(1,2,1)
    imshow(imadjust(C3{iSampleImage})); title('channel #3');
    subplot(1,2,2)
    imshow(BW{iSampleImage});
end

%% calc signal intensity within BW mask
intMean = []; intMedian = []; int95 = []; CFP = []; allCirc = [];

for iPos=1:n
    areas = cat(1,Props{iPos}.Area); %all areas
    circ = cat(1,Props{iPos}.Circularity); %all circularities
    intensityMean = []; intensityMedian = []; intensity95 = []; centroids = []; CFPMedian = []; % place holders
    for iObject=1:length(areas)
        curArea = areas(iObject);
        curCirc = circ(iObject);
        if(curArea>=areaMin & curArea<=areaMax & curCirc >= circCutoff) % test if object is valid by size
            cur = Props{iPos};
            pixelIdxList = cur(iObject).PixelIdxList;
            pixels = C2{iPos}(pixelIdxList); %get the pixels within the mask
            intensityMean = [intensityMean; mean(pixels)]; %get mean YFP within mask
            intensityMedian = [intensityMedian; median(pixels)]; %get median YFP within mask
            intensity95 = [intensity95; prctile(pixels,95)]; %get 95th %ile YFP within mask
            centroids = [centroids; cat(1,cur(iObject).Centroid)]; %find center of mask
            CFPpix = C3{iPos}(pixelIdxList); %get pixels within mask
            CFPMedian = [CFPMedian; median(CFPpix)]; %get median CFP within each mask
        end
        allCirc = [allCirc, curCirc(curArea>=areaMin)]; % save all circularities to plot in a histogram in case the cutoff needs to be adjusted
    
    end
    intMean = [intMean, intensityMean']; 
    intMedian = [intMedian, intensityMedian']; 
    int95 = [int95, intensity95']; 
    CFP = [CFP, CFPMedian'];
end

if tfDebug % plot histogram of cicularities
    figure; hold on;
    histogram(allCirc,'Normalization','probability');
    grid on; box on;
    xlabel('circularity (0-1)')
    ylabel('proportion')
end
%% save variables
fluor = {};
fluor.intMean = intMean;
fluor.intMedian = intMedian;
fluor.int95 = int95;
fluor.CFP = CFP;

end

