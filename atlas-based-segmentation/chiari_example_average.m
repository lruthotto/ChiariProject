%==========================================================================
% 
%  Example for 2D registration of chiari data that averages chiari_example,
%
%   - data                 chiari, omega=(0,1)x(0,1), level=5:8, m=[256,256]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbHyperElastic  
%   - optimization         Gauss-Newton
%
%  This works by summing up the top n closest masks to the reference in
%  terms of SSD and then dividing by n. A threshold is then used to convert
%  this into a binary mask.
% =========================================================================

function vout = chiari_example_average(R, file, varargin)

    if nargin==0
        runMinimalExample; return;
    end
    close all;
    
    %% Initial Setup
    
    % Reference mask along with the threshold and number of samples to
    % average
    Rm         = zeros(256,256);
    thr        = 0.5;
    n          = 10;
    
    % Registration parameters
    omega      = [0,1,0,1];
    m          = [256,256];
    
    % More settings
    %  * Plots refers to to whether it should plot or run silently
    %  * The data set can be changed
    %  * The training size is how much of the data it will actually use
    plots        = 1;
    data         = load('normalizedChiariTrainingData-v2.mat');
    train_size   = 41;
    
    for k=1:2:length(varargin),    % overwrite defaults  
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
    
    % Data
    images       = data.normalTrain(:,:,1:train_size);
    masks        = data.masksTrain(:,:,1:train_size);
    
    % List for storing the SSDs between the reference and all possible
    % templates. The first column counts upwards to keep track of the
    % template id.
    ssd_list  = zeros(train_size, 2);
    ssd_list(:, 1) = 1:train_size;

    %% Calculate SSD for each template image
    
    % Setup reference
    xc = getCellCenteredGrid(omega,m);
    Rc = nnInter(R,omega,xc);

    % Loop through templates and compute SSDs
    for i = 1:train_size
        Ti = images(:,:,i);
        Tc = nnInter(Ti,omega,xc);

        ssd_list(i, 2) = SSD(Tc, Rc, omega, m);
    end
    
    % Sort the list and take the first n
    ssd_sorted = sortrows(ssd_list, 2);
    top_picks = ssd_sorted(1:n,1);
    
    %% Run the image registration and store it in a file for future use
    orient  = @(I) flipud(I)';
    
    % If the file doesn't exist, generate it
    if ~exist(file)
        Tc_list = cell(n, 1);

        for i = 1:n
            % Load the template
            T  = images(:,:,top_picks(i));
            Tm = masks(:,:,top_picks(i));
            
            % Calculate the transformation
            Tc = chiari_example(R, 'T_Tm', {T, Tm}, 'plots', 0);
            
            % Orient and resize it to fit FAIR
            Tc_list{i} = reshape(orient(Tc{1}), [], 1);
        end
        % Save the file for the next run through
        save(file, 'Tc_list')
        
    % If it does exist, load it
    else
        Tc_list = load(file).Tc_list;
        
        % This is in case the file doesn't have enough samples (incomplete
        % file)
        if size(Tc_list) < n
            for i = size(Tc_list)+1:n
                T  = images(:,:,top_picks(i));
                Tm = masks(:,:,top_picks(i));
                
                Tc = chiari_example(R, 'T_Tm', {T, Tm}, 'plots', 0);
                
                Tc_list{i} = reshape(orient(Tc{1}), [], 1);
            end
            save(file, 'Tc_list')
        end
    end
    
    %% Compute transformations and averages
    c_sum = zeros(m(1)*m(2), 1);
    b_sum = zeros(m(1)*m(2), 1);
    
    % Sum up all of the transformed template images
    for i = 1:n
        Tc = Tc_list{i};
        
        c_sum = c_sum + (Tc == 1);
        b_sum = b_sum + (Tc == 2);
    end
    
    % Divide by the total number of images
    c_avg = c_sum / n;
    b_avg = b_sum / n;
    
    % Convert to binary
    c_bin = c_avg > thr;
    b_bin = b_avg > thr;
    
    % Output the mask in the correct format
    Tc = c_bin + 2*b_bin;
    vout{1} = flipud(reshape(Tc, [], 256)');
    
    % Also orient the mask
    Rm = orient(Rm);
    
    %% Plot images
    if plots == 1    
        figure()
        
        subplot(1,3,1)
        viewImage2Dsc(Tc_list{1}, omega, m);
        hold on
        viewContour2D(Rm > 0, omega, m, 'r');
        title("Single T(yc)")
        set(gca,'XTick',[], 'YTick', [])
        axis([0.3    0.9    0.15    0.75]);
        colorbar
        caxis([0 1])
        
        subplot(1,3,2)
        viewImage2Dsc(b_avg + c_avg, omega, m);
        hold on
        viewContour2D(Rm > 0, omega, m, 'r');
        title("Average T(yc)")
        set(gca,'XTick',[], 'YTick', [])
        axis([0.3    0.9    0.15    0.75]);
        colorbar
        caxis([0 1])

        subplot(1,3,3)
        viewImage2Dsc((b_bin + c_bin) > 0, omega, m);
        hold on
        viewContour2D(Rm > 0, omega, m, 'r');
        set(gca,'XTick',[], 'YTick', [])
        title(['Average T(yc) > ' num2str(thr)])
        axis([0.3    0.9    0.15    0.75]);
        colorbar
        caxis([0 1])
    end

    %% Analyze segmentation quality
    if ~isequal(Rm, zeros(256,256))
        Ts = Tc_list{1};
        create_table(Ts, Rm, Tc)
    end
end



%% Function for printing dice
function [] = create_table(Tm, Rm, Tc)
    Tm = reshape(Tm, [], 1);
    Rm = reshape(Rm, [], 1);

    og_b = dice_jaccard(Tm == 1, Rm == 1);
    og_c = dice_jaccard(Tm == 2, Rm == 2);
    og_t = dice_jaccard(Tm > 0,  Rm > 0);

    tn_b = dice_jaccard(Tc == 1, Rm == 1);
    tn_c = dice_jaccard(Tc == 2, Rm == 2);
    tn_t = dice_jaccard(Tc > 0,  Rm > 0);

    T = table([og_b{1}; og_c{1}; og_t{1}], ...
              [tn_b{1}; tn_c{1}; tn_t{1}], ...
              'VariableNames', {'Single Dice', 'Average Dice'}, ...
              'RowName', {'Mask (brain stem)', ...
                          'Mask (cerebellum)', ...
                          'Mask (total)'});

    disp(T);
end


%% Minimal example for chiari_example_average
function runMinimalExample
    id   = 1;
    file = [num2str(id) '_tc.mat'];

    test_data     = load('normalizedChiariTestData-v2.mat');
    reference     = test_data.normalTest(:,:,id);
    mask          = test_data.masksTest(:,:,id);
    
    chiari_example_average(reference, file, 'Rm', mask);
end