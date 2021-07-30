function [biom_bs, biom_cb] = biomarker(mask, dense, varargin)
%% Computes temporal peak spatial average DENSE biomarker
% inputs:
%     mask  = 256x256 matrix of classes of each pixel
%               0: background
%               1: cerebellum
%               2: brainstem
%     dense = 256x256 matrix of peakDisplacement values
% dense_id  = the id number for the DENSE image of interest
% options to include border, see figures, and display calculated biomarkers
%% 
% Establishing Defaults
border = 0;
figs = 0;
display = 1;
for k=1:2:length(varargin)          % overwrite defaults  
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end
if border == 1
    border = blur(dense,1)<2.5;     % defines a border for pixels we will ignore from the average based on too high displacement
    %imagesc(border)
else
    border = ones(256,256);
end

cb = mask == 1;                             % defining cerebellum
dense_cb = dense .* cb .* border;           % finding pixels to average
biom_cb = mean(nonzeros(dense_cb), 'all');  % calculating average

bs = mask==2;                               % defining brain stem
dense_bs = dense .* bs .* border;           % finding pixels to average
biom_bs = mean(nonzeros(dense_bs), 'all');  % calculating average

if figs == 1                %displaying figures of averaged pixels
    figure()
    imagesc(dense_cb)
    
    figure()
    imagesc(dense_bs)
end
if display==1               %displaying calculated biomarkers
    disp(['Cerebellum:' num2str(biom_cb) '   Brain Stem:' num2str(biom_bs)])
end
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
end
output=uint8(output);
end

%% BLUR CODE SOURCE
% https://www.mathworks.com/matlabcentral/answers/450356-how-can-blur-an-image