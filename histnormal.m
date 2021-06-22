function histnormal(refids, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ex:       histnormal([1 2 3 4])   will create 3 figures:
%
%           Figure 1) Original Images
%                     Histograms of Original Images
%                     Cumulative Distributions of Original Images
%           Figure 2) Adjusted Images  (using histeq)
%                     Histograms of Adjusted Images
%                     Cumulative Distributions of Adjusted Images
%           Figure 3) Original Images
%                     Adjusted Images
%           
%           histnormal([1 2 3 4],1) will do the same, and also print the
%           original and adjusted image as original# and adjusted#,
%           respectively, with # corresponding to the image's id number
%
% refids = vector of patient id numbers
%
%  NOTE: this function lines images up side by side in one figure, so doing
%  inputting more than about 5 ids at once (in a vector) is not recommended
%
% Uses chiariTrainingData.mat and MATLAB's Image Processing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
if nargin<2 || isempty(opt)
     opt=0;
end
%%
N = length(refids);
data = load('chiariTrainingData.mat');
images = data.imagesTrain;

%% First figure: 
figure(1)
clf;
for i=1:N
    image = uint8(256 .* (images(:,:,refids(i)) ./ ...  %max adjusted image
        max(max(images(:,:,refids(i))))));
    hold on
    subplot(3,N,i)
    imagesc(image)
    axis equal
    axis off
    colormap gray
    subplot(3,N,i+N)                    %one row below
    histogram(image)
    subplot(3,N,i+2*N)                  %two rows below
    cdfplot(reshape(image,[256*256 1]))
end
max(max(images(:,:,refids(i))))
figure(2)
clf;
%% Second Figure: Equalizing Histogram and Displaying Adjusted Images
for i=1:N
    image = uint8(256 .* (images(:,:,refids(i)) ./ ...
        max(max(images(:,:,refids(i))))));
    hold on
    subplot(3,N,i)
    histeq(image)
    axis equal
    axis off
    colormap gray
    subplot(3,N,i+N)
    histogram(histeq(image))
    subplot(3,N,i+2*N)
    cdfplot(reshape(histeq(image),[256*256 1]))
end
%% Third Figure: Comparing Original and Adjusted Images
figure(3)
clf;
for i=1:N
    image = uint8(256 .* (images(:,:,refids(i)) ./ ...
        max(max(images(:,:,refids(i))))));
    hold on
    subplot(2,N,i)
    imagesc(image)
    axis equal
    axis off
    colormap gray
    subplot(2,N,i+N)
    histeq(image)
end
%% Printing Images as .jpgs
if opt == 1
    for i=1:N
        image = uint8(256 .* (images(:,:,refids(i)) ./ ...     %max adjusted image
            max(max(images(:,:,refids(i))))));
        filename1 = sprintf('%s%d.jpg','original',refids(i));  %filename for original/max-adjusted image
        filename2 = sprintf('%s%d.jpg','adjusted',refids(i));  %filename for adjusted image
        imwrite(image, filename1)
        imwrite(histeq(image), filename2)
    end
end