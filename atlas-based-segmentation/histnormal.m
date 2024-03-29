%==========================================================================
% refids = patient id numbers       vector 
%
% Ex:       histnormal([1 2 3 4], 'figs', 1, 'jpgs', 1, 'mat', 1)  will create 3 figures:
%
%           Figure 1) Original Images
%                     Histograms of Original Images
%                     Cumulative Distributions of Original Images
%           Figure 2) Adjusted Images  (using histeq)
%                     Histograms of Adjusted Images
%                     Cumulative Distributions of Adjusted Images
%           Figure 3) Original Images
%                     Adjusted Images
%
%           It will also print each original and adjusted image to a .jpg
%           file (original#.jpg and adjusted#.jpg) and will store the
%           images in a .mat file (images_original and images_normal, each 
%           256x256xlength(refids) double)
%
%           Defaults with figures on and .jpg/.mat printing off
%
%  NOTE: this function lines images up side by side in one figure, so 
%  inputting more than about 5 ids at once is not recommended for figures
%
% Uses chiariTestingData-v2.mat and MATLAB's Image Processing Toolbox
%==========================================================================

function histnormal(refids, varargin)

    if nargin == 0
         refids=1:51;
    end
    
    %% Defaults
    figs   = 0;
    jpgs   = 0;
    mat    = 0;
    N      = length(refids);
    data   = load('chiariTestData-v2.mat');
    images = data.imagesTest;
    masks  = data.masksTest;
    mean   = data.meanDisplacementTest;
    peak   = data.peakDisplacementTest;

    for k=1:2:length(varargin)    % overwrite defaults  
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end

    %% Figures
    if figs == 1
        % First figure: Original Images and Histograms
        figure(1)       
        clf;
        for i = 1:N
            image = uint8(256 .* (images(:,:,refids(i)) ./ ...  %max adjusted image
                max(max(images(:,:,refids(i))))));
            hold on
            subplot(3,N,i)
            imagesc(image)
            axis equal
            axis off
            colormap gray
            subplot(3,N,i+N)                    %one row below
            histogram(image)
            xlabel('gray value')
            ylabel('number of pixels')
            title('Histogram of Pixel Values')
            subplot(3,N,i+2*N)                  %two rows below
            cdfplot(reshape(image,[256*256 1]))
            axis([0 255 0 1])
        end
        % Second Figure: Equalizing Histogram and Displaying Adjusted Images
        figure(2)
        clf;
        for i = 1:N
            image = uint8(256 .* (images(:,:,refids(i)) ./ ...
                max(max(images(:,:,refids(i))))));
            hold on
            subplot(3,N,i)
            histeq(image)
            axis equal
            axis off
            colormap gray
            subplot(3,N,i+N)
            histogram(histeq(image))
            xlabel('gray value')
            ylabel('number of pixels')
            title('Histogram of Pixel Values')
            subplot(3,N,i+2*N)
            cdfplot(reshape(histeq(image),[256*256 1]))
            axis([0 255 0 1])
        end
        % Third Figure: Comparing Original and Adjusted Images
        figure(3)
        clf;
        for i = 1:N
            image = uint8(256 .* (images(:,:,refids(i)) ./ ...
                max(max(images(:,:,refids(i))))));
            hold on
            subplot(2,N,i)
            imagesc(image)
            axis equal
            axis off
            colormap gray
            colorbar
            title('Original')
            subplot(2,N,i+N)
            histeq(image)
            axis equal
            axis off
            colormap gray
            colorbar
            title('Normalized')
        end
    end
    %% Printing Images as .jpgs
    if jpgs == 1
        for i = refids
            image = uint8(256 .* (images(:,:,refids(i)) ./ ...     %max adjusted image
                max(max(images(:,:,refids(i))))));
            filename1 = sprintf('%s%d.jpg','original',refids(i));  %filename for original/max-adjusted image
            filename2 = sprintf('%s%d.jpg','adjusted',refids(i));  %filename for adjusted image
            imwrite(image, filename1)
            imwrite(histeq(image), filename2)
        end
    end
    %% Creating a .mat file with the new images
    if mat == 1
        originalTest = zeros(256,256,N);
        normalTest = zeros(256,256,N);
        masksTest = zeros(256,256,N);
        meanDisplacementTest = zeros(256,256,N);
        peakDisplacementTest = zeros(256,256,N);
        for i = 1:N
           image = uint8(256 .* (images(:,:,refids(i)) ./ ...     %max adjusted image
                max(max(images(:,:,refids(i))))));
           originalTest(:,:,i) = image;
           normalTest(:,:,i) = histeq(image);
           masksTest(:,:,i) = masks(:,:,refids(i));
           meanDisplacementTest(:,:,i) = mean(:,:,refids(i));
           peakDisplacementTest(:,:,i) = peak(:,:,refids(i));
        end
        clearvars -except originalTest normalTest masksTest meanDisplacementTest peakDisplacementTest
        save('normalizedChiariTestData-v2.mat')
    end

end

