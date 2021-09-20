%==========================================================================
% 
%  Quick script for generating all of the average data. This iterates
%  through all images in the test set
%
% =========================================================================

test_data     = load('normalizedChiariTestData-v2.mat');
    
for i = 1:size(test_data.originalTest, 3)
    R  = test_data.normalTest(:,:,i);
    Rm = test_data.masksTest(:,:,i);

    chiari_example_average(R, [num2str(i) '_tc.mat'], 'Rm', Rm, 'plots', 0);
end