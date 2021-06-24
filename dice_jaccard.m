%==============================================================================
% 
%  Function for calculating the segmentation quality using a soft metric.
%  The matrices must be the same size. They are compared with a combination
%  of intersections and unions.
%  
%  Parameters:
%    - Matrix #1
%    - Matrix #2
%
%  Outputs:
%    - Dice similarity
%    - Jaccard simularity
% ===============================================================================

function varargout = dice_jaccard(varargin)
    if nargin == 0
        runMinimalExample
        return;
    end
    
    data1 = abs(varargin{1});
    data2 = abs(varargin{2});

    % calculate the information
    intersection = sum(data1 .* data2, 'all');
    mask_sum = sum(data1, 'all') + sum(data2, 'all');
    union = mask_sum - intersection;
    
    % edge case if one image is all zeroes
    if mask_sum == 0
        varargout{1} = 1;
        varargout{2} = 1;
        return
    end
        
    % return simularities
    varargout{1} = (2 * intersection) / mask_sum;
    varargout{2} = intersection / union;
    
    return
end

function runMinimalExample
    A = zeros(256, 192);
    B = ones(256, 64);
    data1 = [A B];
    
    A = zeros(256, 128);
    B = ones(256, 128);
    data2 = [A B];
    
    figure()
    colormap(flipud(hot))
    subplot(2, 1, 1)
    imagesc(data1)
    
    subplot(2, 1, 2)
    imagesc(data2)
    
    [d,j] = dice_jaccard(data1, data2);
    disp("dice: " + d);
    disp("jaccard: " + j);
end