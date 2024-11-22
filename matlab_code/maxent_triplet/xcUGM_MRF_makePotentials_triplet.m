function [nodePot,edgePot,tripletPot] = ...
    xcUGM_MRF_makePotentials_triplet...
    (w,nodeMap,edgeMap,tripletMap,tripletStruct)

[nNodes,maxState] = size(nodeMap);
nEdges = tripletStruct.nEdges;
edgeEnds = tripletStruct.edgeEnds;
nStates = tripletStruct.nStates;

nTriplets = tripletStruct.nTriplets;
tripletEnds = tripletStruct.tripletEnds;


nodePot = zeros(nNodes,maxState);
for n = 1:nNodes
    for s = 1:nStates(n)
        if nodeMap(n,s) == 0
            nodePot(n,s) = 1;
        else
            nodePot(n,s) = exp(w(nodeMap(n,s)));
        end
    end
end

if nargout > 1
    edgePot = zeros(maxState,maxState,nEdges);
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        for s1 = 1:nStates(n1)
            for s2 = 1:nStates(n2)
                if edgeMap(s1,s2,e) == 0
                    edgePot(s1,s2,e) = 1;
                else
                    edgePot(s1,s2,e) = exp(w(edgeMap(s1,s2,e)));
                end
            end
        end
    end
end


if nargout > 2
    tripletPot = zeros(maxState,maxState,maxState,nTriplets);
    for t = 1:nTriplets
        n1 = tripletEnds(t,1);
        n2 = tripletEnds(t,2);
        n3 = tripletEnds(t,3);
        for s1 = 1:nStates(n1)
            for s2 = 1:nStates(n2)
                for s3 = 1:nStates(n3)
                    if tripletMap(s1,s2,s3,t) == 0
                        tripletPot(s1,s2,s3,t) = 1;
                    else
                        tripletPot(s1,s2,s3,t) = ...
                            exp(w(tripletMap(s1,s2,s3,t)));
                    end
                end
            end
        end
    end
end








