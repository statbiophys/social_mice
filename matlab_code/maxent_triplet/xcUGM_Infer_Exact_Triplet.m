function  [nodeBel, edgeBel, logZ] = ...
    xcUGM_Infer_Exact_Triplet...
    (nodePot, edgePot, tripletPot, tripletStruct)
% INPUT
% nodePot(node,class)
% edgePot(class,class,edge) where e is referenced by V,E (must be the same
% between feature engine and inference engine)
%
% OUTPUT
% nodeBel(node,class) - marginal beliefs
% edgeBel(class,class,e) - pairwise beliefs
% logZ - negative of free energy

%UGM_assert(prod(double(edgeStruct.nStates)) < 50000000,'Brute Force Exact Inference not recommended for models with > 50 000 000 states');

if tripletStruct.useMex
%      [nodeBel,edgeBel,logZ] = xcUGM_Infer_ExactC(nodePot,edgePot,edgeStruct.edgeEnds,edgeStruct.nStates);

   [nodeBel,edgeBel,logZ] = ...
       xcUGM_Infer_Exact_TripletC ...
       (nodePot,edgePot,tripletStruct.edgeEnds, ...
       tripletPot, tripletStruct.tripletEnds, ...
       tripletStruct.nStates);
else
   [nodeBel,edgeBel,logZ] = Infer_Exact_Triplet...
       (nodePot,edgePot,tripletPot,tripletStruct);
end
end

function  [nodeBel, edgeBel, logZ] = ...
    Infer_Exact_Triplet(nodePot, edgePot, tripletPot, tripletStruct)

[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = tripletStruct.edgeEnds;
tripletEnds = tripletStruct.tripletEnds;
nStates = tripletStruct.nStates;


% Initialize
nodeBel = zeros(size(nodePot));
edgeBel = zeros(size(edgePot));
y = ones(1,nNodes);
Z = 0;
i = 1;
while 1
    
    pot = xcUGM_ConfigurationPotential_triplet...
        (y,nodePot,edgePot,edgeEnds,tripletPot,tripletEnds);
    
    
    % Update Z
    Z = Z + pot;
    
    % Go to next y
    for yInd = 1:nNodes
        y(yInd) = y(yInd) + 1;
        if y(yInd) <= nStates(yInd)
            break;
        else
            y(yInd) = 1;
        end
    end
    
    % Stop when we are done all y combinations
    
      if  yInd == nNodes && y(end) == 1
        break;
    end
end

logZ = log(Z);


end

function assert(pred, str)
% ASSERT Raise an error if the predicate is not true.
% assert(pred, string)

if nargin<2, str = ''; end

if ~pred
  s = sprintf('assertion violated: %s', str);
  error(s);
end
end