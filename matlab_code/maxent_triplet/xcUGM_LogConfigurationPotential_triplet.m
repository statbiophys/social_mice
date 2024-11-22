function [logPot] = xcUGM_LogConfigurationPotential_triplet...
    (y,nodePot,edgePot,edgeEnds,tripletPot,tripletEnds)
% [logPot] = UGM_LogConfigurationPotential(y,nodePot,edgePot,edgeEnds)
nNodes = size(nodePot,1);
nEdges = size(edgeEnds,1);
nTriplets = size(tripletEnds,1);

logPot = 0;

% Nodes
for n = 1:nNodes
   logPot = logPot+log(nodePot(n,y(n)));
   %disp(['node ' num2str(n) ' has logP = ' num2str(log(nodePot(n,y(n))))]);
end
%disp(['logPot after summing up the nodes = ' num2str(logPot)])

% Edges
for e = 1:nEdges
   n1 = edgeEnds(e,1);
   n2 = edgeEnds(e,2);
   logPot = logPot+log(edgePot(y(n1),y(n2),e));
%    if edgePot(y(n1),y(n2),e)~= 1
%        disp([n1, n2, log(edgePot(y(n1),y(n2),e))])
%    end
end

for t = 1:nTriplets
    t1 = tripletEnds(t,1);
    t2 = tripletEnds(t,2);
    t3 = tripletEnds(t,3);
    logPot = logPot + log(tripletPot(y(t1),y(t2),y(t3),t));
end


