%==========================================================================
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
% =========================================================================

function vout = chiari_example(R, varargin)

    if nargin==0
        runMinimalExample; return;
    end
    close all;
    
    %% Initial Setup
    
    % Default values for template, template mask, and reference. Template
    % and template mask are combined so that one can't be inputted without
    % the other
    T_Tm         = {zeros(256,256), zeros(256,256)};
    Rm           = zeros(256,256);
    
    % Registration parameters
    omega        = [0,1,0,1];
    m            = [256,256];
    alpha        = 500;
    mu           = 1;
    lambda       = 0;
    theta        = 1e-2;
    min_level    = 5;
    max_level    = 8;
    
    % More settings
    %  * Plots refers to to whether it should plot or run silently
    %  * The data set can be changed
    %  * The training size is how much of the data it will actually use
    plots        = 1;
    data         = load('normalizedChiariTrainingData-v2.mat');
    train_size   = 41;
    viewPara     = {'viewImage','viewImage2D','colormap','bone(256)'};
    imgPara      = {'imgModel','nnInter'};

    for k=1:2:length(varargin),    % overwrite defaults  
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
    
    % Data
    images       = data.normalTrain(:,:,1:train_size);
    masks        = data.masksTrain(:,:,1:train_size);
    
    % Orient data for FAIR use
    orient = @(I) flipud(I)';
    T_Tm{1} = orient(T_Tm{1});
    T_Tm{2} = orient(T_Tm{2});
    R       = orient(R);
    Rm      = orient(Rm);
    
    % Factor to scale masks when viewing so that they are visible
    mask_scale = 128;
    
    %% Image registration options
    viewImage('reset','viewImage','viewImage2D','colormap', bone(256),'axis','off');
    imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
    trafo('reset','trafo','affine2D');
    regularizer('reset','regularizer','mbHyperElastic','alpha',alpha,'mu',mu,'lambda',lambda);
    
    %% Find the most similar template (SSD)
    if isequal(T_Tm{1}, zeros(256,256)) && isequal(T_Tm{2}, zeros(256,256))
        % Initialize minimum distance and the corresponding index
        min_d = intmax;
        min_i = -1;
        
        % Setup reference
        xc = getCellCenteredGrid(omega,m);
        Rc = nnInter(R,omega,xc);
        
        % Loop through templates, compute SSD, update minimum
        for i = 1:train_size
            Ti = orient(images(:,:,i));
            Tc = nnInter(Ti,omega,xc);

            Dc = SSD(Tc, Rc, omega, m);

            if Dc < min_d
                min_d = Dc;
                min_i = i;
            end
        end
        T_Tm{1} = orient(images(:,:,min_i));
        T_Tm{2} = orient(masks(:,:,min_i));
    end
    
    % More convenient names
    T  = T_Tm{1};
    Tm = T_Tm{2};
    
    %% Get multilevels
    ML = getMultilevel({T,R},omega,m,'fig',2*plots);
    ML_m = getMultilevel({mask_scale*Tm,mask_scale*Rm},omega,m,'fig',2*plots);
    distance('reset','distance','SSD');
    
    %% Calculate the transformation
    yc = MLIR(ML, 'parametric', false,...
        'minLevel', min_level, 'maxLevel', max_level,'plots',plots);
    
    %% Output the transformed mask
    Tc = nnInter(Tm, omega, center(yc, m));
    vout{1} = flipud(reshape(Tc, 256, [])');


    %% Analyze segmentation quality
    if plots, showResults(ML,yc); end
                
    % If a mask is inputted, plot and calculate dice
    if ~isequal(Rm, zeros(256,256))
        % Plotting the transformed mask
        if plots
            showResults(ML_m, yc)

            %Rm and Tm
            FAIRfigure(); clf;
            subplot(1, 2, 1);
            viewImage(R,omega,m);
            hold on
            viewContour2D(Rm, omega, m, 'y');
            viewContour2D(Tm, omega, m, 'b');
            title("R with Rm and Tm")
            axis equal

            %Rm and Tm[y]
            subplot(1, 2, 2);
            viewImage(R,omega,m);
            hold on
            viewContour2D(Rm, omega, m, 'y');
            viewContour2D(Tc, omega, m, 'b');
            title("R with Rm and Tc")
            axis equal
        end
    
        [og_dice, tn_dice] = create_table(Tm, Rm, Tc);
        vout{2} = [og_dice, tn_dice];
    end
end



%% Function for printing dice
function [og_dice, tn_dice] = create_table(Tm, Rm, Tc)
    Tm = reshape(Tm, [], 1);
    Rm = reshape(Rm, [], 1);
    Tc = reshape(Tc, [], 1);

    og_b = dice_jaccard(Tm == 1, Rm == 1);
    og_c = dice_jaccard(Tm == 2, Rm == 2);
    og_t = dice_jaccard(Tm > 0,  Rm > 0);

    tn_b = dice_jaccard(Tc == 1, Rm == 1);
    tn_c = dice_jaccard(Tc == 2, Rm == 2);
    tn_t = dice_jaccard(Tc > 0,  Rm > 0);
    
    og_dice = og_t{1};
    tn_dice = tn_t{1};

    T = table([og_b{1}; og_c{1}; og_t{1}], ...
              [tn_b{1}; tn_c{1}; tn_t{1}], ...
              'VariableNames', {'Original Dice', 'Transformed Dice'}, ...
              'RowName', {'Mask (brain stem)', ...
                          'Mask (cerebellum)', ...
                          'Mask (total)'});

    disp(T);
end



%% FAIR function for viewing the contour of an image
function varargout = viewContour2D(T,omega,m,color,varargin)
    h  = (omega(2:2:end)-omega(1:2:end))./m;
    xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
    [ih, c] = contour(xi(1),xi(2),reshape(T,m)',1, color); axis xy image

    c.LineWidth = 3;

    % the following lines add some nice stuff to the code.
    % if varargin = {'title','FAIR','xlabel','x'}
    % the code evaluates "title('FAIR');xlabel('x');"
    for k=1:2:length(varargin), 
      if ~isempty(varargin{k}), feval(varargin{k},varargin{k+1}); end;
    end;
    if nargout == 1, varargout = {ih};  end;
end



%% Minimal example for chiari_example
function runMinimalExample
    id = 1;

    test_data     = load('normalizedChiariTestData-v2.mat');
    reference     = test_data.normalTest(:,:,id);
    mask          = test_data.masksTest(:,:,id);
    
    chiari_example(reference, 'Rm', mask);
end
%==============================================================================