function [tripleV,tripleE] = xcUGM_makeTripletVE(tripletEnds,nNodes,useMex)
% Calculate : 
% E - [ [indexes of edges connected to node 1][indexes of edges connected to node 2] .... [indexes of edges connected to node nNodes] ]
% V -- V(i) = the sum of the number of edges connected to nodes (1,2,...i-1) plus 1  (the sums of the lengths of the above blocks 1,2..i-1 plus 1)

if nargin < 3
	useMex = 1;
end

if useMex % Note: much more memory efficient
	[tripleV,tripleE] = xcUGM_makeTripletVEC(int32(tripletEnds),int32(nNodes));
else
    %%
	nTriplets = size(tripletEnds,1);
	
	nNei = zeros(nNodes,1);
	nei = zeros(nNodes,0);
	for t = 1:nTriplets
		n1 = tripletEnds(t,1);
		n2 = tripletEnds(t,2);
        n3 = tripletEnds(t,3);
		nNei(n1) = nNei(n1)+1;
		nNei(n2) = nNei(n2)+1;
        nNei(n3) = nNei(n3)+1;
		nei(n1,nNei(n1)) = t;
		nei(n2,nNei(n2)) = t;
        nei(n3,nNei(n3)) = t;
	end
	
	edge_triplet = 1;
	tripleV = zeros(nNodes+1,1,'int32');
	tripleE = zeros(6*nTriplets,1,'int32');
	for n = 1:nNodes
		tripleV(n) = edge_triplet;
		nodeEdges = sort(nei(n,1:nNei(n)));
		tripleE(edge_triplet:edge_triplet+length(nodeEdges)-1,1) = ...
            nodeEdges;
		edge_triplet = edge_triplet+length(nodeEdges);
	end
	tripleV(nNodes+1) = edge_triplet;
end



