%==============================================================================
% 
%  Function for calculating the segmentation quality. The matrices must be
%  the same size. They are compared with a combination of intersections and
%  unions.
%  
%  Parameters:
%    - Matrix #1
%    - Matrix #2
%    - Soft/Hard metric
%
%  Outputs:
%    - Dice similarity
%    - Jaccard simularity
% ===============================================================================

function vout = dice_jaccard(varargin)
    if nargin == 0
        runMinimalExample
        return;
    end
    
    data1 = abs(varargin{1});
    data2 = abs(varargin{2});
    soft = varargin{3};

    smooth = 0.001;
    axes = [1, 2];
        
    % calculate the information
    if(soft == 0)
        mask_sum = sum(data1, axes) + sum(data2, axes);
        inter = sum(data1 .* data2, axes);
        union = mask_sum - inter;
    else
        mask_sum = sum(ceil(data1), axes) + sum(ceil(data2), axes);
        inter = sum(data1 & data2, axes);
        union = sum(data1 | data2, axes);
    end
        
    % return simularities
    vout{1} = mean(2 * (inter + smooth) / (mask_sum + smooth));
    vout{2} = mean((inter + smooth) / (union + smooth));
    
    return
end

function runMinimalExample
    A = zeros(256, 128);
    B = ones(256, 128);
    data1 = [A B];
    
    A = zeros(256, 192);
    B = ones(256, 64);
    data2 = [A B];
    
    figure()
    colormap(flipud(hot))
    subplot(2, 1, 1)
    imagesc(data1)
    
    subplot(2, 1, 2)
    imagesc(data2)
    
    dj = dice_jaccard(data1, data2, 1);
    disp("Dice: " + dj(1));
    disp("Jacc: " + dj(2));
end