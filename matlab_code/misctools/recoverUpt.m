function [output] = recoverUpt(input, matrixSize)
% Given the upper triangular part of a matrix in the format of 1d array,
% recover the matrix

    output = zeros(matrixSize);

    counter = 1;
    for i = 1:matrixSize
        for j = (i+1):matrixSize
            output(i,j) = input(counter);
            counter = counter + 1;
        end
    end
    
    output = output + output';

end
