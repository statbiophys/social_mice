function [cc3, c3] = compute_cc3 (s, p)
%% Compute triplet

%%
[nNodes, nSamples] = size(s);

[~, m, ~, c2] = compute_cc (s, s, p, p);

% Compute correlation tensor
% if esFlag % equal state correlation
%     c3 = zeros(nNodes, nNodes, nNodes);
%     for i = 1:nNodes
%         for j = 1:nNodes
%             for k = 1:nNodes
%                     c3(i,j,k) = ...
%                         mean((s(i,:)==s(j,:)).*(s(j,:)==s(k,:)));
%             end
%         end
%     end
% else
%     
    
c3 = zeros(nNodes, nNodes, nNodes, p, p, p);
cc3 = zeros(nNodes, nNodes, nNodes, p, p, p);



for i = 1:nNodes
    for j = 1:nNodes
        for k = 1:nNodes
            for ii = 1:p
                for jj = 1:p
                    for kk = 1:p
                        c3(i,j,k,ii,jj,kk) = ...
                            mean((s(i,:)==ii).*...
                            (s(j,:)==jj).*...
                            (s(k,:)==kk));
                        cc3(i,j,k,ii,jj,kk) = ...
                            c3(i,j,k,ii,jj,kk)...
                            -m(i,ii)*c2(j,k,jj,kk) ...
                            -m(j,jj)*c2(i,k,ii,kk) ...
                            -m(k,kk)*c2(j,i,jj,ii) ...
                            +2*m(i,ii)*m(j,jj)*m(k,kk);
                    end
                end
            end
        end
    end
end


end

