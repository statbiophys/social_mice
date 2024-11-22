function [samples] = ...
    xcUGM_Sample_Gibbs_triplet...
    (nodePot,edgePot,tripletPot,tripletStruct,burnIn,y)
% [samples] = UGM_Sample_Gibbs(nodePot,edgePot,edgeStruct,burnIn,y)
% Single Site Gibbs Sampling

if nargin < 6
% Initialize
%%
[junk y] = max(nodePot,[],2);
end

if tripletStruct.useMex
    samples = ...
        xcUGM_Sample_Gibbs_tripletC(nodePot,edgePot,tripletPot,...
        tripletStruct.edgeEnds, tripletStruct.nStates, tripletStruct.V,...
        tripletStruct.E, ...
        tripletStruct.tripletEnds, tripletStruct.tripleV,...
        tripletStruct.tripleE, ...
        int32(tripletStruct.maxIter), int32(burnIn), int32(y));
else
    samples = ...
        Sample_Gibbs_triplet(nodePot,edgePot,tripletPot,tripletStruct,burnIn,y);
end

end

function [samples] = Sample_Gibbs_triplet...
    (nodePot,edgePot,tripletPot,tripletStruct,burnIn,y)

%%
[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = tripletStruct.edgeEnds;
V = tripletStruct.V; % careful, here is where V and E enter
E = tripletStruct.E;

tripletEnds = tripletStruct.tripletEnds;
tripleV = tripletStruct.tripleV; % careful, here is where V and E enter
tripleE = tripletStruct.tripleE;
nTriplets = tripletStruct.nTriplets;

nStates = tripletStruct.nStates;
maxIter = tripletStruct.maxIter;

samples = zeros(nNodes,0);
%%
for i = 1:burnIn+maxIter
    for n = 1:nNodes
    %%
        % Compute Node Potential
        pot = nodePot(n,1:nStates(n));

        % Find Neighbors
        edges = E(V(n):V(n+1)-1);

        % Multiply Edge Potentials
        for e = edges(:)'
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);

            if n == edgeEnds(e,1)
                ep = edgePot(1:nStates(n1),y(n2),e)';
            else
                ep = edgePot(y(n1),1:nStates(n2),e);
            end
            pot = pot .* ep;
        end

        % Find triple edges
        triplet_edges = tripleE(tripleV(n):tripleV(n+1)-1);
        for t = triplet_edges(:)'
            n1 = tripletEnds(t,1);
            n2 = tripletEnds(t,2);
            n3 = tripletEnds(t,3);

            % pay attention to dimension. maybe need to squeeze then
            % transpose
            if n == n1 %edgeEnds(e,1)
                tp = tripletPot(1:nStates(n1),y(n2),y(n3),t)';
            elseif n == n2
                tp = tripletPot(y(n1),1:nStates(n2),y(n3),t)';
%                 tp = edgePot(y(n1),1:nStates(n2),e);
            else
                tp = squeeze(tripletPot(y(n1),y(n2),1:nStates(n3),t))';
            end
            pot = pot .* tp;
        end

        % Sample State;
        y(n) = sampleDiscrete(pot./sum(pot));
    end
    
    if i > burnIn
        samples(:,i-burnIn) = y;
    end
end
end