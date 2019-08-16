% Copyright (c) 2019 Surjyendu Ray
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%% Function to generate frequency of mutations matrix from mutation counts data
function [F_from_SampleData, scaling_factor_float] =  transform_elkebir_input_data_into_F_matrix(input_file)

    fileF = input_file;
    fID = fopen(fileF,'r');

    delimiterIn = '\t';
    headerlinesIn = 1;
    input_mat = importdata( fileF,delimiterIn,headerlinesIn );
    disp( input_mat );

    data_table = input_mat.data;
    [n,T] = size(data_table);
    numberNodes = n;
    samples = T/2;
    
    F_from_SampleData = [];
    
    %Left set shall be the counts of the reference allele
    %right set shall be the counts of the alternate allele
    left_set = data_table(:,1:2:end);
    right_set = data_table(:,2:2:end);

    first_half = nan(1,samples);
    
    for t = 1:numberNodes
        for i = 1:samples
             
             %so F_in_reality = 2 *  right_counts / (right_counts + left_counts)

             % we threshold the values at 1 because we know that they
             % cannot be larger than 1.
             mut_freq_GPU = (2 * right_set(t,i)) / (left_set(t,i) + right_set(t,i));
             if (isnan(mut_freq_GPU)) % if there is a nan in left or right. we assume f will be zero
                mut_freq_GPU = 0;
             end
             mut_freq_GPU = min(1,mut_freq_GPU);
             
             first_half(i) = mut_freq_GPU;
        end   
       
        %Returning the F_matrix represented by the first_half array
        F_from_SampleData = [F_from_SampleData; first_half];
    end   
    
      
    %Getting a scaling factor by summing every pair of columns and then
    %averaging
    sum_col_pairs = data_table(:,1:2:end) + data_table(:,2:2:end);
    scaling_factor_float = mean(sum_col_pairs(~isnan(sum_col_pairs(:))));
        
    fprintf('The scaling factor for Ancestree simulation = %d\n', scaling_factor_float);

    fclose(fID);
end