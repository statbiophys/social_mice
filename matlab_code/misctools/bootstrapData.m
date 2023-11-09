function [output, outSample] = bootstrapData (input, binSize, prcKeep, tau)
% XC 20170329
% [output, outSample] = bootstrapData (input, binSize, prcKeep)
% prcKeep : precentage of the data that I want to use. Default: 0.5
% binSize: default = 50
% tau (added 20181008): in dynamical model, the time difference

if nargin < 4
    tau = 0;
end

if nargin < 3
    prcKeep = 0.5;
end

if nargin < 2
%     disp('Use binsize = 50.');
    binSize = 50;
end

[nSample,~] = size(input);
if nSample == 1
    input = input';
    [nSample,~] = size(input);
end

% rng('shuffle'); 


kk = randperm(floor(nSample/binSize));
% kk = kk(1:(floor(nSample/(binSize*2))-1));
kk = kk(1:floor(nSample*prcKeep/binSize));
index = [];
for j = 1:length(kk)
%     index = [index (kk(j)-1)*binSize+(1:binSize)];
    index = [index (kk(j)-1)*binSize+(1:binSize-tau)];
end

output = input(index,:);

otherIdx = setdiff(1:nSample, index);
outSample = input(otherIdx,:);

end
