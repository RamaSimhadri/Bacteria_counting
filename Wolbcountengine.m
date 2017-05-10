% This script will count the number of Wolbachia bacteria in a confocal Z-stack. 
% FIJI or ImageJ have scripts for edge detection and counting objects, 
% but this script does better when the bacteria are clumped together.
% Once the clumps are identified, local maxima, i.e. the centers of each 
% bacterium within the clump are counted, which produces a more accurate
% estimate of the number of bacteria in the image.
%
% Input: Confocal Z-stacks of cell lines
% Output: Plots of the bacteria identified and the number of bacteria in each plane

close all
clear all
clc

% Read through the 1st to 7th slice in the Z-stack
S = 1;
N = 7;
Idata = cell(N,1);
loadcnt = 1;
for i = S:N
    
    Idata{loadcnt} = imread(strcat('1_C001Z00',num2str(i),'.tif'));
    loadcnt = loadcnt+1;
    
end

map = hot(16);
pixelthreshold = 300; % VERY IMPORTANT, this determines the removal of noise
normmax = 16; %
suborder = 1; % How many images from top and from below the stackt to subtract
intvec = [2 8 14];

scrsz = get(0,'ScreenSize');
figimg = figure('Position',[100 100 scrsz(3)-200 scrsz(4)-200]);

bs = S+1:N-1; %slices to analyze
Lf = cell(length(bs),1);
dim1 = cell(length(bs),1);
dim2 = cell(length(bs),1);
Istore = cell(length(bs),1);
Mplots = length(bs);
Nplots = 2;
cnt = 1;
for k = bs

    if k == 1
        disp('Cannot analyze. Need lower layer');
        close(figimg)
        break;
    elseif k == N
        disp('Cannot analyze. Need upper layer');
        close(figimg)
        break;
    else
        remv = [k-suborder:k-1 k+1:k+suborder];
    end
    
    base = k;
    I = Idata{base};%load image 

    for j = remv
        I = I - Idata{j};
    end

    %filters to remove noise from the background
    Ltmp = medfilt2(I,[6 6]);
    G2 = fspecial('average',[5 5]);
    W2 = imfilter(Ltmp,G2,'replicate');
    H = fspecial('gaussian',[3 3], 0.8);
    L2 = imfilter(W2,H,'replicate');
    L2db = double(L2);
    L2 = L2db/max(max(L2db))*normmax;
    
    Lf{cnt} = L2;
    
    dim1{cnt} = size(Lf{cnt},1);
    dim2{cnt} = size(Lf{cnt},2);
    
    figure(figimg)
    subplot(Nplots,Mplots,Mplots+cnt)
    imshow(L2,map)
    xlabel('Processed Image')
    
    cnt = cnt+1;
    
end

Lindividual = cell(length(bs),1);
NP = cell(length(bs),1);
for k = 1:length(bs)
    
    % identifying edges and counting the number of "bacterial" clumps
    K = bwconncomp(Lf{k});
    numPixels = cellfun(@numel,K.PixelIdxList);
    bacteria = find(numPixels > pixelthreshold);
    NP{k} = numPixels(bacteria);
    lumps = cell(length(bacteria),1);
    lumpcnt = 1;
    Lindividual{k} = zeros(dim1{k},dim2{k});
    for indx = 1:length(bacteria)
        lumps{lumpcnt} = K.PixelIdxList{bacteria(indx)};
        Lindividual{k}(lumps{lumpcnt}) = Lf{k}(lumps{lumpcnt});
        lumpcnt = lumpcnt+1;
    end
    
    subplot(Nplots,Mplots,k)
    imshow(Lindividual{k},map)
    
    xlabel(['Image number: ',num2str(bs(k))])
    hold on
    
    L3 = Lindividual{k};
        
    lp = 1;
    [B, L] = bwboundaries(L3> intvec(lp),'noholes');
    for kl = 1:length(B)
        plot(B{kl}(:,2),B{kl}(:,1),'.','Color','b')
    end
    Bf = B;
    nele = length(B);
    disp(nele)
    % identifying individual bacteria within the clumps 
    for kl = 1:length(Bf)
        
        Bf2 = Bf;
        for lp = 2:length(intvec)
            ll = intvec(lp);
            bndclr = [ll/normmax 1-ll/normmax ll/normmax];
            xlow = min(Bf2{kl}(:,1));
            xhigh = max(Bf2{kl}(:,1));
            ylow = min(Bf2{kl}(:,2));
            yhigh = max(Bf2{kl}(:,2));
            L4 = L3(xlow:xhigh,ylow:yhigh);
            [Bf3, L] = bwboundaries(L4 > intvec(lp));

            for kl1 = 1:length(Bf3)
                plot(ylow+Bf3{kl1}(:,2),xlow+Bf3{kl1}(:,1),'.','Color',bndclr)
            end

            if length(Bf3) > 1
                nele = nele - 1 + length(Bf3);
            end

        end

    end
    
    title(['Number of Lumps: ',num2str(length(B)),' \newline Number of Bacteria: ',num2str(nele)])
end

