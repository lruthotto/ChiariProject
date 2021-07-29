function [biom_bs, biom_cb] = biomarker(mask, dense_id, opt) %mask will eventually be an input
%% 
% inputs:
%     mask  = matrix of classes of each pixel
%               0: background
%               1: cerebellum
%               2: brainstem
% dense_id  = the id number for the DENSE image of interest
% opt for border (1) or no border (0)
% If you would like to 
%% 
if nargin<3 || isempty(opt)
     opt = 0;
end
data=load('chiariTrainingData-v2.mat');
dense=data.peakDisplacementTrain(:,:,dense_id);
%mask=data.masksTrain(:,:,dense_id);   %% Temporary -- eventually would be inputting a mask

if opt == 1
    border = blur(dense,1)<2.5;     %uses dense data to make a hard border for where we see CSF
    imagesc(border)
else
    border = ones(256,256);
end

% figure(1)
% clf;
cb = mask==1;
dense_cb = dense .* cb .* border;
% imagesc(dense_cb)
biom_cb = mean(nonzeros(dense_cb), 'all');

% figure(2)
% clf;
bs = mask==2;    
dense_bs = dense .* bs .* border;
% imagesc(dense_bs)
biom_bs = mean(nonzeros(dense_bs), 'all');

disp(['Cerebellum:' num2str(biom_cb) '   Brain Stem:' num2str(biom_bs)])
end

function [output] = blur(A,w)
[row col] = size(A);
A=uint8(A);
B=nan(size(A) + (2*w));
B(w+1:end-w,w+1:end-w)=A;
output = 0*A;
for i=w+1:row+w
  for j=w+1:col+w
    tmp=B(i-w:i+w,j-w:j+w);
    output(i-w,j-w)=mean(tmp(~isnan(tmp)));
  end
output=uint8(output);
end

%% BLUR CODE SOURCE
% https://www.mathworks.com/matlabcentral/answers/450356-how-can-blur-an-image