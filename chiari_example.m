%==============================================================================
% 
%  Example for 2D registration of chiari data,
%
%   - data                 hands, omega=(0,20)x(0,25), level=3:7, m=[128,128]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mfElastic
%   - optimization         Gauss-Newton
% ===============================================================================


%% Initial Setup
close all

% Parameters
Template_ID = 18;
Reference_ID = 28;

omega     = [0,20,0,25];
m         = [256,256];
viewPara = {'viewImage','viewImage2D','colormap','bone(256)'};

% Data
data = load('chiariTrainingData.mat');

images = data.imagesTrain;
masks = data.masksTrain;

orient = @(I) flipud(I)';
mask_scale = 128;

dataT = orient(images(:,:,Template_ID));
dataR = orient(images(:,:,Reference_ID));
dataT_mask = mask_scale .* orient(masks(:,:,Template_ID));
dataR_mask = mask_scale .* orient(masks(:,:,Reference_ID));

% Multilevel representation
ML = getMultilevel({dataT,dataR},omega,m,'fig',2);
ML_mask = getMultilevel({dataT_mask,dataR_mask},omega,m);

% More options
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',1e3,'mu',1,'lambda',0);

%% Build a weighted SSD distance that ignores background
xc = reshape(getCellCenteredGrid(omega,m),[],2);
c  = (omega(2:2:end)+omega(1:2:end))./2;
dataW = exp(- .01*(xc(:,1)-c(1)).^2 - .01*(xc(:,2)-c(2)).^2)+.1;

%% Plot the reference and the weights
figure(1);clf;
subplot(1,2,1);
viewImage(dataR,omega,m);
title('dataR - reference image');
colorbar

subplot(1,2,2);
ph = viewImage2Dsc(dataW,omega,m,'caxis',[0 1]);
title('dataW - weights');
colorbar;

%% Generate multilevel representation of data images and weights
[ML,minLevel,maxLevel] = getMultilevel({dataT,dataR,dataW},omega,m,'names',{'T','R','W'});

% extract multilevel representation of weights
%
% 1) get discrete data of weighs from ML
% 2) get coefficients and evaluate image model on grids
% 3) store weights for each level as diagonal matrix 
% 4) provide multi-level representation of weights to distance
%

MLw = cell(maxLevel,1);
for lvl=minLevel:maxLevel
    WD = ML{lvl}.W;  % 1) 
    WC = WD;         % 2) coefficients for linear inter
    Wc = linearInter(WC,omega,getCellCenteredGrid(omega,ML{lvl}.m)); % 3)
    MLw{lvl}.Wc = diag(sparse(Wc)); 
    MLw{lvl}.m  = ML{lvl}.m;
    ML{lvl} = rmfield(ML{lvl},'W');
end

distance('reset','distance','SSD','weights',MLw);

%% Calculate and display the transformation

[yc,wc,his] = MLIR(ML,'parametric',1,'plotIter',1,'plotMLiter',1);

% Also apply the transformation to the mask
showResults(ML,yc)
showResults(ML_mask, yc)

%==============================================================================