%% Multiple cell analysis
% scans over multiple intensity thresholds to identify cells in each image
% in a sequence
% Can be used to identify a sequence with an individual cell or with
% multiple cells over time
inputScale = inputdlg('Enter scale val in pixels per micron');
scaleVal = str2double(inputScale{1});
startThresh = 280; % dimmest cell thresh value
threshValList = logspace(log10(startThresh),log10(10000),100); % analysis scans through these possible thresholds to identify cells
[~, imPath] = uigetfile('*.tiff', 'Choose first file of image sequence for fluorescent cells');
sequence = dir(fullfile(imPath,'*.tiff'));
if isempty(sequence)
    sequence = dir(fullfile(imPath,'*.tif'));
end
imFiles = {sequence.name};
for i = 1:length(imFiles)
    imFiles{i} = fullfile(imPath, imFiles{i});
end
numIms = length(imFiles);
% identify and track individual cells
figure
colormap gray
cellAnalysis = struct('loc',[],'MFI',[]); % for storing data by cell
numCells = 0;
BGVals = zeros(1,numIms);
for i = 1:numIms
    im = double(imread(imFiles{i}));
    imagesc(im)
    imForBG = im;
    if i == 1
        BGTest = modalIntensity(imForBG);
    end
    imForBG(im>1.1*BGTest) = 0;
    BGVals(i) = modalIntensity(imForBG);
    hold on
    for j = 1:length(threshValList)
        binIm = im > threshValList(j); % cur cell thresh
        binIm = imfill(binIm, 'holes');
        connComps = bwconncomp(binIm);
        areaProps = regionprops(connComps,'Area','Centroid','Circularity','PixelIdxList');
        pxIdxList = connComps.PixelIdxList;
        areas = [areaProps.Area];
        circVals = [areaProps.Circularity];
        targetArea = 40*scaleVal^2;
        keepAreas = areas > 0.5*targetArea & areas < 2*targetArea;% & circVals > 0.5;
        %     keepAreas = areas > 500 & areas < 5000;
        areaProps = areaProps(keepAreas);
        if j == 1
            allAreaProps = areaProps;
        else
            for m = 1:length(areaProps)
                placed = false;
                for n = 1:length(allAreaProps)
                    curLoc = areaProps(m).Centroid;
                    compLoc = allAreaProps(n).Centroid;
                    curArea = areaProps(m).Area;
                    compArea = allAreaProps(n).Area;
                    if norm(curLoc-compLoc,2) < 2.5*scaleVal
                        placed = true;
                        if abs(curArea-targetArea) < abs(compArea-targetArea)
                            allAreaProps(n) = areaProps(m);
                            break;
                        end
                    end
                end
                if ~placed
                    allAreaProps(length(allAreaProps)+1) = areaProps(m);
                end
            end
        end
    end
    areaProps = allAreaProps;
    if isempty(areaProps)
        continue
    else
        for j = 1:length(areaProps)
            curLoc = areaProps(j).Centroid;
            curPx = areaProps(j).PixelIdxList;
            emptyIm = zeros(512,512);
            emptyIm(curPx) = 1;
            curBound = bwboundaries(emptyIm);
            %             selPx = curRatioIm(pxIdxList{j});
            % select pixels from a radius of 5 microns
            [xMat,yMat] = meshgrid(1:512,1:512);
            selLogic = sqrt((xMat-curLoc(1)).^2 + (yMat-curLoc(2)).^2) <= 5*scaleVal;
            selPx = im(selLogic);
            orderedPx = sort(selPx,'ascend');
            
            numPx = round(40*scaleVal^2); % adjust this depending on the size of the cell
            if length(selPx) < numPx+6
                curMFI = mean(orderedPx(1:end-5));
            else
                curMFI = mean(orderedPx(end-numPx-5+1:end-5));
            end
            placed = false;
            if numCells > 0 && i > 1
                for k = 1:length(cellAnalysis)
                    prevLoc = cellAnalysis(k).loc(i-1,:);
                    centroidInCurBody = inpolygon(prevLoc(1),prevLoc(2),curBound{1}(:,2),curBound{1}(:,1));
                    if centroidInCurBody%norm(curLoc - prevLoc) < 2.5*scaleVal  % close -> consider same cell
                        cellAnalysis(k).loc(i,:) = curLoc;
                        cellAnalysis(k).MFI(i) = curMFI;
                        placed = true;
                        break
                    end
                end
            end
            if ~placed % then not identified with cell initiated prev
                numCells = numCells + 1;
                cellAnalysis(numCells).loc = nan(numIms,2);
                cellAnalysis(numCells).MFI = nan(numIms,1);
                cellAnalysis(numCells).loc(i,:) = curLoc;
                cellAnalysis(numCells).MFI(i,:) = curMFI;
            end
        end
    end
    for k = 1:length(areaProps)
        curCentroid = areaProps(k).Centroid;
        plot(curCentroid(1),curCentroid(2),'r*')
    end
    title(sprintf('Cells at frame %d',i))
    drawnow
    hold off
end
BG = mean(BGVals);

%%
keepLogic = true(numCells,1);
for i = 1:numCells
    if sum(~isnan(cellAnalysis(i).MFI)) < 15
        keepLogic(i) = false;
    end
end
cellAnalysis = cellAnalysis(keepLogic);

%% make plots
figure
hold on
for i = 1:length(cellAnalysis)
    curMFI = cellAnalysis(i).MFI;
%     curMFI = curMFI(~isnan(curMFI));
        plot((curMFI-BG))
end