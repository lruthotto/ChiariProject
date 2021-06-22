function histnormal(refids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays the max-adjusted image, histogram, and cumulative distribution 
% for inputted reference id numbers in figure 1, with the normalized image,
% histogram, and cumulative distribution in figure 2.  Figure 3 prints
% the original and adjusted images next to each other, for user comparison.
%
% refids = vector of patient id numbers
% Uses chiariTrainingData.mat and MATLAB's Image Processing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
N = length(refids);
data = load('chiariTrainingData.mat');
images = data.imagesTrain;

%%
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
