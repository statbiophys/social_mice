function [my_estimate] = xc_UGM_histogram_method_triplet...
    (Y,tripletStruct,tripletMap,edgeMap,nodeMap,...
    param0,param1)

% Goal: Apply histogram method to estimate the expectation value of a
% function
% This is specific to UGM, i.e. undirected graphs
%
%
% Xiaowen Chen
% 2018-10-30

%%
nInstances = size(Y,1);

dparam = param1 - param0;

nParams = length(dparam);

my_denominator = 0.0;
my_numerator = zeros(nParams, 1);

for i = 1:nInstances
    y = Y(i,:);
    
    a = xcUGM_MRF_computeSuffStat_triplet...
        (y,nodeMap,edgeMap,tripletMap,tripletStruct);
    my_exponent = dot(dparam, a);
    my_denominator = my_denominator + exp(my_exponent);
    my_numerator = my_numerator + exp(my_exponent)*a;
end

my_estimate = my_numerator./my_denominator;


%%
% [nInstances, nNodes] = size(Y);
% 
% dparam = param1 - param0;
% 
% nParams = length(dparam);
% 
% my_denominator = 0.0;
% my_numerator = zeros(nParams, 1);
% 
% for i = 1:nInstances
%     y = Y(i,:);
%     
%     a = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);
%     my_exponent = dot(dparam, a);
%     my_denominator = my_denominator + exp(my_exponent);
%     my_numerator = my_numerator + exp(my_exponent)*a;
% end
% 
% my_estimate = my_numerator./my_denominator;
% 
% if edgeStruct.useMex
% %     my_estimate = xc_UGM_histogram_methodC(edgeStruct.edgeEnds, edgeStruct.nStates,...
% %         edgeStruct.V, edgeStruct.E, edgeMap, nodeMap, int32(Y),);
%     my_estimate = xc_histogram_method(Y,edgeStruct,edgeMap,nodeMap,param0,param1);
% else
%     my_estimate = xc_histogram_method(Y,edgeStruct,edgeMap,nodeMap,param0,param1);
% end



end
% 
% function [my_estimate] = xc_histogram_method(Y,edgeStruct,edgeMap,nodeMap,param0,param1)
% 
% nInstances = size(Y,1);
% 
% dparam = param1 - param0;
% 
% nParams = length(dparam);
% 
% my_denominator = 0.0;
% my_numerator = zeros(nParams, 1);
% 
% for i = 1:nInstances
%     y = Y(i,:);
%     
%     a = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);
%     my_exponent = dot(dparam, a);
%     my_denominator = my_denominator + exp(my_exponent);
%     my_numerator = my_numerator + exp(my_exponent)*a;
% end
% 
% my_estimate = my_numerator./my_denominator;
% 
% end
