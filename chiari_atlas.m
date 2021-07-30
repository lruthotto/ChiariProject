%==========================================================================
% 
%  Quick script that makes it convenient to access chiari_example and
%  chiari_example_average
%
% =========================================================================

function chiari_atlas(R_ID, type)
    training_size = 0;
    test_data     = load('normalizedChiariTestData-v2.mat');
    
    R  = test_data.normalTest(:,:,training_size + R_ID);
    Rm = test_data.masksTest(:,:,training_size + R_ID);
    
    if isequal(type, 'sgl')
        chiari_example(R, 'Rm', Rm)
    elseif isequal(type, 'avg')
        chiari_example_average(R, [num2str(R_ID) '_Tc.mat'], 'Rm', Rm)
    end
