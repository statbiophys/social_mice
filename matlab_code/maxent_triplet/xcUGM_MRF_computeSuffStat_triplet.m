function [suffStat] = ...
    xcUGM_MRF_computeSuffStat_triplet...
    (Y,nodeMap,edgeMap,tripletMap,tripletStruct)

%%
[nNodes,maxState] = size(nodeMap);
nEdges = tripletStruct.nEdges;
edgeEnds = tripletStruct.edgeEnds;
% nStates = tripletStruct.nStates;

nTriplets = tripletStruct.nTriplets;
tripletEnds = tripletStruct.tripletEnds;

nParams = max([nodeMap(:);edgeMap(:);tripletMap(:)]);

nInstances = size(Y,1);
suffStat = zeros(nParams,1);

for i = 1:nInstances
   y = Y(i,:);

   %%
   for n = 1:nNodes
      if nodeMap(n,y(n)) > 0
         suffStat(nodeMap(n,y(n))) = ...
             suffStat(nodeMap(n,y(n))) + 1;
      end
   end
   %%

   for e = 1:nEdges
      n1 = edgeEnds(e,1);
      n2 = edgeEnds(e,2);
      if edgeMap(y(n1),y(n2),e) > 0
         suffStat(edgeMap(y(n1),y(n2),e)) = ...
             suffStat(edgeMap(y(n1),y(n2),e)) + 1;
      end
   end

   for t = 1:nTriplets
       n1 = tripletEnds(t,1);
       n2 = tripletEnds(t,2);
       n3 = tripletEnds(t,3);
       if tripletMap(y(n1),y(n2),y(n3),t) > 0
           suffStat(tripletMap(y(n1),y(n2),y(n3),t)) = ...
               suffStat(tripletMap(y(n1),y(n2),y(n3),t)) + 1;
       end
   end    


end




