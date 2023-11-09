function [output, index_new, rset] = random_cyclic_shuffle(input, given_index)

% if size(input,1) > size(input,2)
%     input = input';
% end
[nNodes, nSamples] = size(input);

rset = [];

if nargin < 2
    output = input;
    index_old = repmat(1:nSamples,nNodes,1);
    index_new = index_old;
    
    for i = 1:nNodes
        r = randi(nSamples);
        rset = [rset r];
        output(i,r+1:end) = input(i,1:end-r);
        output(i,1:r) = input(i,end-r+1:end);

        index_new(i,r+1:end) = index_old(i,1:end-r);
        index_new(i,1:r) = index_old(i,end-r+1:end);
    end
else
    for i = 1:nNodes
        output(i,:) = input(i,given_index(i,:));
    end
    index_new = given_index;
end

end