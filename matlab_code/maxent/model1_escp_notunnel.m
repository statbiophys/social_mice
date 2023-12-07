function [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates)
% Define the model that encode pairwise interaction in the same box and
% that each mouse has its own box preference
% escp = equal state, constant parameters


maxState = nStates;
maxIter = 1e4;

%% construct the graphic model using UGM

adj = ones(nNodes,nNodes)-diag(ones(nNodes,1));
edgeStruct = UGM_makeEdgeStruct(adj, nStates);
edgeStruct.maxIter=maxIter;
nodeMap=zeros(nNodes,maxState,'int32');
for iState=1:nStates-1
    nodeMap(:,iState)=1+nNodes*(iState-1):nNodes*iState;
end
nEdges=edgeStruct.nEdges;
edgeMap=zeros(maxState,maxState,nEdges,'int32');
for e = 1:edgeStruct.nEdges
    for iq = 1:nStates
        edgeMap(iq,iq,e)=e+nNodes*(nStates-1);
    end
end



end