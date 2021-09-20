%==========================================================================
%  Function for calculating the segmentation quality. The matrices must be
%  the same size. They are compared with a combination of intersections and
%  unions.
%  
%  Parameters:
%    - Matrix #1
%    - Matrix #2
%
%  Outputs:
%    - Dice similarity
%    - Jaccard simularity
% =========================================================================

function vout = dice_jaccard(data1, data2, varargin)

    if nargin == 0
        runMinimalExample
        return;
    end

    smooth = 0.001;
    
    for k=1:2:length(varargin),    % overwrite defaults  
        eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;
        
    % calculate the information
    mask_sum = sum(ceil(data1), 'all') + sum(ceil(data2), 'all');
    inter = sum(data1 & data2, 'all');
    union = sum(data1 | data2, 'all');
        
    % return simularities
    vout{1} = mean(2 * (inter + smooth) / (mask_sum + smooth));
    vout{2} = mean((inter + smooth) / (union + smooth));
    
    return
end



function runMinimalExample
    A = ones(256, 128);
    B = zeros(256, 128);
    data1 = [A B];
    
    A = ones(256, 64);
    B = zeros(256, 192);
    data2 = [A B];
    
    figure()
    colormap(flipud(hot))
    subplot(2, 1, 1)
    imagesc(data1)
    
    subplot(2, 1, 2)
    imagesc(data2)
    
    dj = dice_jaccard(data1, data2);
    disp("Dice: " + dj(1));
    disp("Jacc: " + dj(2));
end