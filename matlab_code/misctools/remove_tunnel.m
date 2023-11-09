function [s_notunnel] = remove_tunnel(s)
% combine the tunnel state with the state that follows

nNodes = size(s,1);
% s = squeeze(s2(1,:,:));

tmax = [];
for imice = 1:nNodes
    mys = s(imice,:);
    tmax = [tmax find(mys>0, 1, 'last')];
end

% tmax = min(tmax)

% s = s(:,1:tmax);
% nT = tmax;

%%
s_notunnel = s;
for imice = 1:nNodes
    
    mys = s(imice,:);
    mys_new = mys;
    for t = tmax(imice):-1:2
        if mys_new(t-1) == 0
            mys_new(t-1) = mys_new(t);
        end
    end
    s_notunnel(imice,:) = mys_new;
end

s_notunnel = s_notunnel(:,1:min(tmax));

%%

end