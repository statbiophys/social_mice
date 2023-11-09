% XC 20170306
% Given the upper triangular part of a matrix in the format of 1d array,
% recover the matrix

function [output] = makeUpt(input)
%     
%     output = triu(input,1);
%     output = reshape(output', length(input)^2,1);
%     output = output(output>0);
    
    l=size(input,1);
    output = zeros(l*(l-1)/2,1);
    count = 1;
    for i = 1:l-1
        for j = i+1:l
            output(count) = input(i,j);
            count = count+1;
        end
    end

end