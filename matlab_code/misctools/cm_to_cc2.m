function [cc2] = cm_to_cc2(cij, mir, nNodes, nStates)
% 2020-01-19
% Utility code, computing connected two point cu
%%
if numel(cij) == nNodes*(nNodes-1)/2
    if min(size(cij))==1
        cij = recoverUpt(cij, nNodes);
    end

    if min(size(mir))==1
        mir = reshape(mir, nNodes, nStates-1);
    end

    mir = cat(2, mir, 1-sum(mir,2));

    cc2 = zeros(nNodes);
    cc2 = cij - mir*mir';

    cc2 = makeUpt(cc2);
elseif numel(cij) == nNodes*(nNodes-1)*(nStates-1)*(nStates-1)/2
    cij = reshape(cij, nStates-1, nStates-1, nNodes*(nNodes-1)/2);
    
    mir = reshape(mir, nNodes, nStates-1);
    mir = cat(2, mir, 1-sum(mir,2));

    cc2 = zeros(size(cij));
    count = 1;
    for i = 1:nNodes
        for j = i+1:nNodes
            cc2(:,:,count) = cij(:,:,count) - mir(i,1:4)'*mir(j,1:4);
            
            count = count + 1;
        end
    end
else
    cc2 = nan(nNodes);
end



end