function [result, strArray] = MNI_valGrab_v2(inFile, roidir, fxnHandles, csvname,patID,scantype,procfile,functionsInput)
% takes a volume and computes statistics for the voxels in various regions

% input:
    % inFile (char) : path to an image volume in 2mm MNI space
    % roidir (char) : path to a directory containing right-sided structural masks
        % MNI_2mm_[region]_R.nii.gz (or just .nii)
    % fxnHandles (cell array of function handles) : the functions you want
        % to apply to the voxels in regions. {@mean @std} is typical
        % note that a cell array is specified with curly braces
    % csvname (char) : filename to which results will be written
    
% returns:
    % strArray: (m+1) x n string array, m = # of fxn handles, n = # of masks in roidir * 2
        % row 1: region name
        % row n+1: the value of the fxnHandles{n} applied to that region
    % result: a double array equivalent to strArray but without the first (header) row
    % also writes the result to csvname
    % note that region names are sorted alphabetically, and are presented
    % in pairs (left and right)

%% need to add the path to load_nii
% addpath(genpath('/Users/manusdonahue/Documents/Sky/MNI_valGrab/nifti_tools')); % this should be the folder that load_nii is in
    
%% check what roi files we have and sort them
roidirRegex = strcat(roidir,'/MNI_2mm_*_R.nii*'); % can use NiFTi or zipped NiFTi
listing = dir(roidirRegex);

files = [string('')];
regions = [string('')];
for i = 1:length(listing)
    currentFile = string(getfield(listing(i),'name'));
    files(i) = strcat(roidir, '/', currentFile);
    
    split = strsplit(currentFile, '_');
    currentRegion = split(3); % all files should be structured as MNI_2mm_[region]_R.nii.gz
    regions(i) = currentRegion;
end

[sorted_regions, sorted_indices] = sort(regions);
sorted_files = files(sorted_indices);

%% load our volume of interest
[volume] = load_nii(inFile,[],[],[],[],[],1);
volumeMatrix = double(volume.img);

%% go through each mask and calculate statistics

outputRegionNames = [string('')];
outputStats = zeros(length(fxnHandles), length(sorted_regions)*2);
for i = 1:length(sorted_regions)
    workingFile = sorted_files(i);
    regionBasename = sorted_regions(i);
    rName = strcat(regionBasename, '_R');
    lName = strcat(regionBasename, '_L');
    
    outputRegionNames(i*2-1) = lName;
    outputRegionNames(i*2) = rName;
    
    rightMask = load_nii(char(workingFile));
    rightMaskMat = double(rightMask.img);
    leftMaskMat = flipud(rightMaskMat);
    
    maskedRight = volumeMatrix;
    maskedRight(rightMaskMat == 0) = NaN;
    maskedRight(maskedRight == 0) = NaN;
    
    rightClean = rmmissing(reshape(maskedRight, prod(size(maskedRight)), 1)); % this flattens the array and removes nans, which are outside the mask
    
    maskedLeft = volumeMatrix;
    maskedLeft(leftMaskMat == 0) = NaN;
    maskedLeft(maskedLeft == 0) = NaN;
    leftClean = rmmissing(reshape(maskedLeft, prod(size(maskedLeft)), 1)); % this flattens the array and removes nans, which are outside the 
    
%     for l=1:size(maskedLeft,3)
%         imagesc(volumeMatrix(:,:,l))
%         title(workingFile)
%         pause
%         imagesc(leftMaskMat(:,:,l))
%         title(workingFile)
%         pause
%     end
    
    % iterate through the functions specified and store the results
    for j = 1:length(fxnHandles)
        
        fxn = fxnHandles{j};
        %disp('----')
        %disp('Calculating')
        %disp(fxn)
        %disp('for')
        %disp(regionBasename)
        %disp('----')
        
        rightResult = fxn(rightClean);
        leftResult = fxn(leftClean);
        
        outputStats(j, i*2-1) = leftResult;
        outputStats(j, i*2) = rightResult;
    end
end

%% write results

strArray = [(outputRegionNames); string(outputStats)];
result = outputStats;

dataFile=fopen(csvname,'a');

s = dir(csvname);
if s.bytes == 0
    fwrite(dataFile,['patID, scan type, processed file,functionsInput,', char(strjoin(outputRegionNames, ', ')), newline]);
end

theSize = size(outputStats);
for i = 1:theSize(1)
    disp(i)
    disp(outputStats(i,:))
    outputStats(i,(isnan(outputStats(i,:))))=-999;
    fwrite(dataFile, [patID,',',scantype,',',procfile,',',functionsInput,',',char(strjoin(string(outputStats(i,:)), ', ')), newline]);
end

fclose(dataFile);

%% testing
% inFile = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/examples/rOEF.nii.gz';
% roidir = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/MNI_valGrab_masks';
% fxnHandles = {@mean @std};
% csvname = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/rOEF_valgrab_results.csv';
% 
% result = MNI_valGrab(inFile, roidir, fxnHandles, csvname);
