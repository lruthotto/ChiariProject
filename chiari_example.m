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
function[] = chiari_example(Template_ID, Reference_ID) 

%% Defaults 
if nargin<2 || isempty(Reference_ID) 
     Reference_ID = 28; 
end 
if nargin<1 || isempty(Template_ID) 
     Template_ID = 18; 
end 

%% Initial Setup
close all
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
distance('reset','distance','SSD');
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mbHyperElastic','alpha',1e3,'mu',1,'lambda',0);

%% Calculate and display the transformation
yc = MLIR(ML, 'parametric', false,...
    'minLevel', 5, 'maxLevel', 8,'plots',1);

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
ycc = center(yc, m);

mg_yc = nnInter(dataT, omega, ycc);
mg_r = reshape(dataR, [], 1);
mk_yc = round(nnInter(dataT_mask, omega, ycc) ./ 128);
mk_r = round(reshape(dataR_mask, [], 1) ./ 128);

mg_t = reshape(dataT, [], 1);
mk_t = round(reshape(dataT_mask, [], 1) ./ 128);

max_value = max([mg_yc; mg_r]);

og_mg_data = {mg_t ./ max_value, mg_r ./ max_value};
og_mk_c_data = {mk_t == 1, mk_r == 1};
og_mk_b_data = {mk_t == 2, mk_r == 2};
og_mk_t_data = {mk_t & mk_t, mk_r & mk_r};

tn_mg_data = {mg_yc ./ max_value, mg_r ./ max_value};
tn_mk_c_data = {mk_yc == 1, mk_r == 1};
tn_mk_b_data = {mk_yc == 2, mk_r == 2};
tn_mk_t_data = {mk_yc & mk_yc, mk_r & mk_r};

og_mg_s = dice_jaccard(og_mg_data{1}, og_mg_data{2}, 0);
og_mk_b = dice_jaccard(og_mk_b_data{1}, og_mk_b_data{2}, 1);
og_mk_c = dice_jaccard(og_mk_c_data{1}, og_mk_c_data{2}, 1);
og_mk_t = dice_jaccard(og_mk_t_data{1}, og_mk_t_data{2}, 1);

tn_mg_s = dice_jaccard(tn_mg_data{1}, tn_mg_data{2}, 0);
tn_mk_b = dice_jaccard(tn_mk_b_data{1}, tn_mk_b_data{2}, 1);
tn_mk_c = dice_jaccard(tn_mk_c_data{1}, tn_mk_c_data{2}, 1);
tn_mk_t = dice_jaccard(tn_mk_t_data{1}, tn_mk_t_data{2}, 1);

T = table([og_mg_s{1}; og_mk_b{1}; og_mk_c{1}; og_mk_t{1}], ...
          [tn_mg_s{1}; tn_mk_b{1}; tn_mk_c{1}; tn_mk_t{1}], ...
          [og_mg_s{2}; og_mk_b{2}; og_mk_c{2}; og_mk_t{2}], ...
          [tn_mg_s{2}; tn_mk_b{2}; tn_mk_c{2}; tn_mk_t{2}], ...
          'VariableNames',{'Original Dice', 'Transformed Dice', 'Original Jaccard', 'Transformed Jaccard'}, ...
          'RowName', {'Magnitude', ...
                      'Mask (brain stem)', ...
                      'Mask (cerebellum)', ...
                      'Mask (total)'});

disp(T);
% table2latex(T, 'table');
%==============================================================================