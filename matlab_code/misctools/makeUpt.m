function [output] = makeUpt(input)
% Given the upper triangular part of a matrix in the format of 1d array,
% recover the matrix

    
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
