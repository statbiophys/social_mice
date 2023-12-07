function [output] = histogramEqualization(input)

    output=input;
    nsample=length(input);
    [~,sortIdx]=sort(input);
    
    sortedOutput=linspace(0,1,nsample);
    
    for i = 1:nsample
        output(sortIdx(i))=sortedOutput(i);
    end

end