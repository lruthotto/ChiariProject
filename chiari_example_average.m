function [] = chiari_example_average(Reference_ID)
    %% Initial Setup
    close all
    if nargin < 1, Reference_ID = 1; end
    
    % Data
    data = load('normalizedChiariTraining.mat');
    images = data.images_normal;
    masks = data.masks;
    
    omega     = [0,1,0,1];
    m         = [256,256];
    
    avg_w_rg  = 1;
    avg_thr   = 0.6;
    n         = 30;
    d_size    = size(images, 3);
    
    ssd_list = zeros(size(images, 3), 2);
    ssd_list(:, 2) = 1:d_size;
    
    orient = @(I) flipud(I)';
    dataR = orient(images(:,:,Reference_ID));
    dataR_mask = orient(masks(:,:,Reference_ID));

    %% Calculate SSD for each template image
    for i = 1:d_size
        dataT = orient(images(:,:,i));

        xc = getCellCenteredGrid(omega,m);
        Tc = nnInter(dataT,omega,xc);
        Rc = nnInter(dataR,omega,xc);

        ssd_list(i, 1) = SSD(Tc, Rc, omega, m);
    end
    
    ssd_sorted = sortrows(ssd_list, 1);
    top_picks = ssd_sorted(2:n+1,2);
    
    %% Run the image registration and store it in a file for future use
    if ~exist([num2str(Reference_ID) '_yc.mat'])
        yc_list = cell(n, 1);

        for i = 1:n
            c = top_picks(i);
            yc_list(i) = chiari_example(Reference_ID, c, 'plots', 0);
        end
        save([num2str(Reference_ID) '_yc.mat'], 'yc_list')
    else
        yc_list = load([num2str(Reference_ID) '_yc.mat']).yc_list;
        if size(yc_list) < n
            for i = size(yc_list):n
                c = top_picks(i);
                yc_list(i) = chiari_example(Reference_ID, 'Template_ID', c, 'plots', 0);
            end
            save([num2str(Reference_ID) '_yc.mat'], 'yc_list')
        end
    end
    
    %% Compute transformations and averages
    avg_weights = logspace(0 + avg_w_rg, 0 - avg_w_rg, n);
    
    c_sum = zeros(m(1)*m(2), 1);
    b_sum = zeros(m(1)*m(2), 1);
    for i = 1:n
        yc = yc_list{i};
        dataT_mask = orient(masks(:,:,top_picks(i)));
        
        Tc = nnInter(dataT_mask, omega, center(yc, m));
        
        c_sum = c_sum + avg_weights(i) .* (Tc == 1);
        b_sum = b_sum + avg_weights(i) .* (Tc == 2);
    end

    cb_max = max(max(c_sum), max(b_sum));
    
    c_avg = c_sum ./ cb_max;
    b_avg = b_sum ./ cb_max;
    
    c_bin = c_avg > avg_thr;
    b_bin = b_avg > avg_thr;
    
    %% plot images
    figure()
    
    subplot(2,2,1)
    viewImage2Dsc(256 * b_avg, omega, m);
    hold on
    viewContour2D(dataR_mask == 2, omega, m);
    axis([0.3    0.9    0.15    0.75]);
    colorbar
    
    subplot(2,2,2)
    viewImage2Dsc(256 * c_avg, omega, m);
    hold on
    viewContour2D(dataR_mask == 1, omega, m);
    axis([0.3    0.9    0.15    0.75]);
    colorbar
    
    subplot(2,2,3)
    viewImage2Dsc(256 * b_bin, omega, m);
    hold on
    viewContour2D(dataR_mask == 2, omega, m);
    axis([0.3    0.9    0.15    0.75]);
    colorbar
    
    subplot(2,2,4)
    viewImage2Dsc(256 * c_bin, omega, m);
    hold on
    viewContour2D(dataR_mask == 1, omega, m);
    axis([0.3    0.9    0.15    0.75]);
    colorbar
end