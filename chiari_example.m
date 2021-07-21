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
function vout = chiari_example(Reference_ID, varargin) 
    %% Initial Setup
    close all;
    omega        = [0,1,0,1];
    m            = [256,256];
    default_ref  = 6;
    
    viewPara     = {'viewImage','viewImage2D','colormap','bone(256)'};
    imgPara      = {'imgModel','nnInter'};
    
    Template_ID  = -1;
    plots        = 1;
    weighted     = 0;
    weight_scale = 25;
    min_level    = 5;
    max_level    = 8;

    for k=1:2:length(varargin),    % overwrite defaults  
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
    
    % Data
    data = load('normalizedChiariTraining.mat');
    images = data.images_normal;
    masks = data.images_masks;
    
    orient = @(I) flipud(I)';
    mask_scale = 128;

    if nargin < 1, Reference_ID = default_ref; end
    dataR = orient(images(:,:,Reference_ID));
    dataR_mask = orient(masks(:,:,Reference_ID));
    
    %% Image registration options
    viewImage('reset','viewImage','viewImage2D','colormap', bone(256),'axis','off');
    imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
    trafo('reset','trafo','affine2D');
    regularizer('reset','regularizer','mbHyperElastic','alpha',500,'mu',1,'lambda',0);
    
    %% Calculate weights
    if weighted
        FAIRfigure([]); clf;
        viewImage(dataR, omega, m)
        focus_center = ginput(1);
        close;
        
        xc = reshape(getCellCenteredGrid(omega,m),[],2);
        c  = focus_center;
        w  = weight_scale;
        dataW = exp(- w*(xc(:,1)-c(1)).^2 - w*(xc(:,2)-c(2)).^2)+.1;
        
        Wt = nnInter(dataW,omega,getCellCenteredGrid(omega,m));
        Wc = diag(sparse(Wt));
    end
    
    %% Find the most similar template
    if Template_ID == -1       
        min_dist = intmax;
        min_index = -1;
        
        for i = 1:size(images, 3)
            if i == Reference_ID
                continue
            end
            
            dataT = orient(images(:,:,i));
    
            xc = getCellCenteredGrid(omega,m);
            Tc = nnInter(dataT,omega,xc);
            Rc = nnInter(dataR,omega,xc);
            
            if weighted
                Dc = SSD(Tc, Rc, omega, m, 'weights', Wc);
            else
                Dc = SSD(Tc, Rc, omega, m);
            end

            if Dc < min_dist
                min_dist = Dc;
                min_index = i;
            end
        end
        
        Template_ID = min_index;
    end

    dataT = orient(images(:,:,Template_ID));
    dataT_mask = orient(masks(:,:,Template_ID));
    
    %% Get multilevels
    if weighted == 0
        ML = getMultilevel({dataT,dataR},omega,m,'fig',2*plots);
        ML_mask = getMultilevel({mask_scale .* dataT_mask, mask_scale .* dataR_mask},omega,m,'fig',2*plots);
        distance('reset','distance','SSD');
    else
        %% Plot the weights
        figure(); clf;
        subplot(1,2,1);
        viewImage(dataR,omega,m);
        title('dataR - reference image');
        colorbar

        subplot(1,2,2);
        viewImage2Dsc(dataW,omega,m,'caxis',[0 1])
        title('dataW - weights');
        colorbar;

        %% Generate multilevel representation of data images and weights
        ML = getMultilevel({dataT,dataR,dataW},omega,m,'names',{'T','R','W'},'fig',2*plots);
        ML_mask = getMultilevel({mask_scale .* dataT_mask, mask_scale .* dataR_mask},omega,m,'fig',2*plots);

        % Extract multilevel representation of weights
        %
        % 1) get discrete data of weighs from ML
        % 2) get coefficients and evaluate image model on grids
        % 3) store weights for each level as diagonal matrix 
        % 4) provide multi-level representation of weights to distance

        MLw = cell(max_level,1);
        for lvl=min_level:max_level
            WD = ML{lvl}.W;  % 1) 
            WC = WD;         % 2) coefficients for linear inter
            Wc = nnInter(WC,omega,getCellCenteredGrid(omega,ML{lvl}.m)); % 3)
            MLw{lvl}.Wc = diag(sparse(Wc)); 
            MLw{lvl}.m  = ML{lvl}.m;
            ML{lvl} = rmfield(ML{lvl},'W');
        end

        distance('reset','distance','SSD','weights',MLw);
    end
    
    %% Calculate and display the transformation
    yc = MLIR(ML, 'parametric', false,...
        'minLevel', 5, 'maxLevel', 8,'plots',plots);
    
    Tc = nnInter(dataT_mask, omega, center(yc, m));
    
    vout{1} = Tc;
    
    % Also apply the transformation to the mask
    if plots
        showResults(ML,yc)
        showResults(ML_mask, yc)
    end
    
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

%     %Use the transformation function on the template mask
%     ycc = center(yc, m); %change to a centered grid
%     [Tmask] = nnInter(dataT_mask,omega,ycc);
%     
%     figure(); clf;
%     %Reference image
%     viewImage(dataR,omega,m);
%     hold on
%     %Plot transformed mask
%     viewContour2D(Tmask, omega, m);
%     axis equal
%     colorbar

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