function [tripletStruct] = ...
    xcUGM_makeTripletStruct(adj, adj_triplet,nStates,useMex,maxIter)
% [edgeStruct] = UGM_makeEdgeStruct(adj,nStates,useMex,maxIter)
%
% adj - nNodes by nNodes adjacency matrix (0 along diagonal)
%
% Combine edgeStruct into tripletStruct

if nargin < 4
    useMex = 1;
end
if nargin < 5
    maxIter = 100;
end

%%
nNodes = int32(length(adj_triplet));
[i j] = ind2sub([nNodes nNodes],find(adj));
nEdges = length(i)/2;
edgeEnds = zeros(nEdges,2,'int32');
eNum = 0;
for e = 1:nEdges*2
   if j(e) < i(e)
       edgeEnds(eNum+1,:) = [j(e) i(e)];
       eNum = eNum+1;
   end
end

[V,E] = UGM_makeEdgeVE(edgeEnds,nNodes,useMex);


tripletStruct.edgeEnds = edgeEnds;
tripletStruct.V = V; % vertex % I need to figure out what those means, V and E
tripletStruct.E = E; % edges
tripletStruct.nEdges = size(edgeEnds,1);


%%
[i j k] = ind2sub([nNodes nNodes nNodes],find(adj_triplet));
% nEdges = length(i)/2;
nTriplets = length(i)/6; % check this
tripletEnds = zeros(nTriplets,3,'int32');
% edgeEnds = zeros(nEdges,2,'int32');
eNum = 0;

% for e = 1:nEdges*2
%    if j(e) < i(e)
%        edgeEnds(eNum+1,:) = [j(e) i(e)];
%        eNum = eNum+1;
%    end
% end

for t = 1:nTriplets*6
   if j(t) < i(t)
       if k(t) < j(t)
           tripletEnds(eNum+1,:) = [k(t) j(t) i(t)];
           eNum = eNum+1;
       end
   end
end


[tripleV,tripleE] = xcUGM_makeTripletVE(tripletEnds,nNodes,useMex); % what does this mean?


tripletStruct.tripletEnds = tripletEnds;
tripletStruct.tripleV = tripleV; % vertex % I need to figure out what those means, V and E
tripletStruct.tripleE = tripleE; % edges
tripletStruct.nNodes = nNodes;
tripletStruct.nTriplets = size(tripletEnds,1);

% Handle other arguments
if isscalar(nStates)
   nStates = repmat(nStates,[nNodes 1]);
end
tripletStruct.nStates = int32(nStates(:));
tripletStruct.useMex = useMex;
tripletStruct.maxIter = int32(maxIter);


% edgeStruct.edgeEnds = edgeEnds;
% edgeStruct.V = V;
% edgeStruct.E = E;
% edgeStruct.nNodes = nNodes;
% edgeStruct.nEdges = size(edgeEnds,1);
% 
% % Handle other arguments
% if isscalar(nStates)
%    nStates = repmat(nStates,[nNodes 1]);
% end
% edgeStruct.nStates = int32(nStates(:));
% edgeStruct.useMex = useMex;
% edgeStruct.maxIter = int32(maxIter);


