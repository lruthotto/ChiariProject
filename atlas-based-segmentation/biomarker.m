function [biom_bs, biom_cb] = biomarker(mask, dense, varargin)
%% Computes temporal peak spatial average DENSE biomarker
% inputs:
%     mask  = 256x256 matrix of classes of each pixel
%               0: background
%               1: cerebellum
%               2: brainstem
%     dense = 256x256 matrix of peakDisplacement values
% dense_id  = the id number for the DENSE image of interest
% options to see figures and display calculated biomarkers
%% 
% Establishing Defaults
figs = 0;
display = 1;
for k=1:2:length(varargin)          % overwrite defaults  
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

cb = mask == 1;                             % defining cerebellum
dense_cb = dense .* cb;                     % finding pixels to average
biom_cb = mean(nonzeros(dense_cb), 'all');  % calculating average

bs = mask==2;                               % defining brain stem
dense_bs = dense .* bs;                     % finding pixels to average
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