%==============================================================================
% 
%  Example for 2D registration of chiari data,
%
%   - data                 chiari, omega=(0,20)x(0,25), level=3:7, m=[256,256]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     rigid2D
%   - regularizer          mbElastic
%   - optimization         Gauss-Newton
% ===============================================================================


%% Initial Setup
close all

% Parameters
Template_ID = 18;
Reference_ID = 28;

omega     = [0,20,0,25];
m         = [256,256];
viewPara  = {'viewImage','viewImage2D','colormap','bone(256)'};
imgPara   = {'imgModel','linearInter'};

% Data
data = load('chiariTrainingData.mat');
normData = load('normalizedChiariTraining.mat');

images = normData.images_normal;
masks = data.masksTrain;

orient = @(I) flipud(I)';
mask_scale = 128;

dataT = orient(images(:,:,Template_ID));
dataR = orient(images(:,:,Reference_ID));
dataT_mask = mask_scale .* orient(masks(:,:,Template_ID));
dataR_mask = mask_scale .* orient(masks(:,:,Reference_ID));

% Multilevel representation
viewImage('reset',viewPara{:});
imgModel('reset',imgPara{:});

ML = getMultilevel({dataT,dataR},omega,m,'fig',2);
ML_mask = getMultilevel({dataT_mask,dataR_mask},omega,m);

% More options
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','NGF');
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mbElastic','alpha',1e3,'mu',1,'lambda',0);

%% Calculate and display the transformation
[yc,wc,his] = MLIR(ML, 'parametric', false,...
    'minLevel', 3, 'maxLevel', 8,'plots',1);

% Also apply the transformation to the mask
showResults(ML,yc)
showResults(ML_mask, yc)

%% Plotting the transformed mask

% %Use the transformation function on the template mask
% ycc = center(yc, m); %change to a centered grid
% [Tmask] = nnInter(dataT_mask,omega,ycc);
% 
% figure(); clf;
% %Reference image
% viewImage(dataR,omega,m);
% hold on
% %Plot transformed mask
% viewContour2D(Tmask, omega, m);
% axis equal
% colorbar

%% Soft metric of segmentation quality
[T_mg, R_mg] = imgModel('coefficients',ML{length(ML)}.T,ML{length(ML)}.R,omega,'out',0);
[T_mk, R_mk] = imgModel('coefficients',ML_mask{length(ML)}.T,ML_mask{length(ML)}.R,omega,'out',0);
xc = getCellCenteredGrid(omega, m);
ycc = center(yc, m);

mg_yc = imgModel(T_mg, omega, ycc);
mg_r = imgModel(R_mg, omega, xc);
mk_yc = imgModel(T_mk, omega, ycc);
mk_r = imgModel(R_mk, omega, xc);

max_value = max([mg_yc; mg_r]);
model_mg_yc = mg_yc ./ max_value;
model_mg_r = mg_r ./ max_value;

model_mk_yc = mk_yc & mk_yc;
model_mk_r = mk_r & mk_r;

[d_mg,j_mg] = dice_jaccard(model_mg_yc, model_mg_r, 0);
[d_mk,j_mk] = dice_jaccard(model_mk_yc, model_mk_r, 1);

disp("dice (mag):  " + d_mg);
disp("jacc (mag):  " + j_mg);
disp("dice (mask): " + d_mk);
disp("jacc (mask): " + j_mk);

%==============================================================================