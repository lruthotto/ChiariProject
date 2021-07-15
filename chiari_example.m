%==============================================================================
% 
%  Example for 2D registration of chiari data,
%
%   - data                 chiari, omega=(0,1)x(0,1), level=5:8, m=[256,256]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbHyperElastic  
%   - optimization         Gauss-Newton
% ===============================================================================
function [] = chiari_example(Reference_ID, Template_ID) 
    %% Initial Setup
    close all
    omega     = [0,1,0,1];
    m         = [256,256];
    viewPara  = {'viewImage','viewImage2D','colormap','bone(256)'};
    imgPara   = {'imgModel','linearInter'};

    % Data
    data = load('normalizedChiariTraining.mat');
    images = data.images_normal;
    masks = data.masks;
    
    orient = @(I) flipud(I)';
    mask_scale = 128;
    
    %% Find the most similar template
    if nargin < 2
        if nargin < 1, Reference_ID = 28; end
        
        dataR = orient(images(:,:,Reference_ID));
        
        min_dist = intmax;
        min_index = -1;
        
        for i = 1:52
            if i == Reference_ID
                continue
            end
            
            dataT = orient(images(:,:,i));
    
            xc = getCellCenteredGrid(omega,m);
            Tc = nnInter(dataT,omega,xc);
            Rc = nnInter(dataR,omega,xc);

            Dc = SSD(Tc, Rc, omega, m);
            
            if Dc < min_dist
                min_dist = Dc;
                min_index = i;
            end
        end
        disp(min_dist)
        Template_ID = min_index;
    end

    dataT = orient(images(:,:,Template_ID));
    dataR = orient(images(:,:,Reference_ID));
    dataT_mask = orient(masks(:,:,Template_ID));
    dataR_mask = orient(masks(:,:,Reference_ID));
    
    %% Image registration options
    viewImage('reset','viewImage','viewImage2D','colormap', bone(256),'axis','off');
    imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
    distance('reset','distance','SSD');
    trafo('reset','trafo','affine2D');
    regularizer('reset','regularizer','mbHyperElastic','alpha',500,'mu',1,'lambda',0);

    %% Multilevel registration
    ML = getMultilevel({dataT,dataR},omega,m);
    ML_mask = getMultilevel({mask_scale .* dataT_mask, mask_scale .* dataR_mask},omega,m);

    %% Calculate and display the transformation
    yc = MLIR(ML, 'parametric', false,...
        'minLevel', 5, 'maxLevel', 8,'plots',1);
    
    % Also apply the transformation to the mask
    showResults(ML,yc)
    showResults(ML_mask, yc)
    
    disp("Reference: " + Reference_ID)
    disp("Template:  " + Template_ID)

    %% Compute NGF
%     xc = getCellCenteredGrid(omega,m);
%     Tc = nnInter(min(dataT_mask, 1), omega, xc);
%     Tc_y = nnInter(min(dataT_mask, 1), omega, center(yc, m));
%     Rc = nnInter(dataR,omega,xc);
% 
%     Dc = NGF(Tc, Rc, omega, m, 'edge', 1e-6);
%     Dc_y = NGF(Tc_y, Rc, omega, m, 'edge', 1e-6);
%
%     disp("NGF Before: " + Dc)
%     disp("NGF After:  " + Dc_y)
    
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

    %% Analyze segmentation quality
    ycc = center(yc, m);

    create_table(dataT_mask, dataR_mask, omega, ycc)
    
end

%% Function for printing dice/jaccard
function [] = create_table(dataT_mask, dataR_mask, omega, ycc)
    mk_yc = nnInter(dataT_mask, omega, ycc);
    mk_t = reshape(dataT_mask, [], 1);
    mk_r = reshape(dataR_mask, [], 1);

    og_mk_c_data = {mk_t == 1, mk_r == 1};
    og_mk_b_data = {mk_t == 2, mk_r == 2};
    og_mk_t_data = {mk_t & mk_t, mk_r & mk_r};

    tn_mk_c_data = {mk_yc == 1, mk_r == 1};
    tn_mk_b_data = {mk_yc == 2, mk_r == 2};
    tn_mk_t_data = {mk_yc & mk_yc, mk_r & mk_r};

    og_mk_b = dice_jaccard(og_mk_b_data{1}, og_mk_b_data{2}, 1);
    og_mk_c = dice_jaccard(og_mk_c_data{1}, og_mk_c_data{2}, 1);
    og_mk_t = dice_jaccard(og_mk_t_data{1}, og_mk_t_data{2}, 1);

    tn_mk_b = dice_jaccard(tn_mk_b_data{1}, tn_mk_b_data{2}, 1);
    tn_mk_c = dice_jaccard(tn_mk_c_data{1}, tn_mk_c_data{2}, 1);
    tn_mk_t = dice_jaccard(tn_mk_t_data{1}, tn_mk_t_data{2}, 1);

    T = table([og_mk_b{1}; og_mk_c{1}; og_mk_t{1}], ...
              [tn_mk_b{1}; tn_mk_c{1}; tn_mk_t{1}], ...
              [og_mk_b{2}; og_mk_c{2}; og_mk_t{2}], ...
              [tn_mk_b{2}; tn_mk_c{2}; tn_mk_t{2}], ...
              'VariableNames',{'Original Dice', 'Transformed Dice', 'Original Jaccard', 'Transformed Jaccard'}, ...
              'RowName', {'Mask (brain stem)', ...
                          'Mask (cerebellum)', ...
                          'Mask (total)'});

    disp(T);
end
%==============================================================================