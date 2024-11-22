function [tripletStruct, nodeMap, edgeMap, tripletMap] = ...
    model5_escp_notunnel_triplet(nNodes, nStates)
% xiaowen chen 20200-1-13
% escp = equal state, constant parameters


maxState = nStates;
maxIter = 1e4;

%% construct the graphic model using UGM

adj = ones(nNodes,nNodes)-diag(ones(nNodes,1));


%% Now let's try making a map + struct for triplet interaction

adj_triplet = ones(nNodes,nNodes,nNodes);
for i = 1:nNodes
    for j = 1:nNodes
        for k = 1:nNodes
            if k == i || i == j || j == k
                adj_triplet(i,j,k)=0;
            end
        end
    end
end

%%
useMex = 0;
tripletStruct = xcUGM_makeTripletStruct(adj, adj_triplet, nStates, useMex);
tripletStruct.maxIter = maxIter;

%%
% edgeStruct = UGM_makeEdgeStruct(adj, nStates);
% edgeStruct.maxIter = maxIter;
nodeMap = zeros(nNodes,maxState,'int32');
for iState = 1:nStates-1
    nodeMap(:,iState) = 1+nNodes*(iState-1):nNodes*iState;
end

nEdges = tripletStruct.nEdges;
edgeMap = zeros(maxState,maxState,nEdges,'int32');
for e = 1:tripletStruct.nEdges
    for iq = 1:nStates
        edgeMap(iq,iq,e) = e+nNodes*(nStates-1);
    end
end


% nodeMap=zeros(nNodes,maxState,'int32');
% for iState=1:nStates-1
%     nodeMap(:,iState)=1+nNodes*(iState-1):nNodes*iState;
% end
nTriplets = tripletStruct.nTriplets;
tripletMap = zeros(maxState,maxState,maxState,nTriplets,'int32');
for t = 1:tripletStruct.nTriplets
    for iq = 1:nStates
        tripletMap(iq,iq,iq,t) = ...
            t+nNodes*(nStates-1)+tripletStruct.nEdges;
    end
end





end