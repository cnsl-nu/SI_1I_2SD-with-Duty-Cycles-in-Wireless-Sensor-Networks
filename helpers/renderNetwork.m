function fig = renderNetwork(Net, fig, tit, x0)


nodes = Net(6);
nodes = nodes{1};
U = triu(Net{5}{1});
[ii,jj] = find(U);
xx = nodes(1,:);
yy = nodes(2,:);
fig = figure(fig);

for k=1:length(ii)
    xxx = [xx(ii(k)),xx(jj(k))];
    yyy = [yy(ii(k)),yy(jj(k))];
    hold on
    p1 = plot(xxx,yyy,'-k');
    p1.Color(4) = 0.5;
end
% s = scatter(xx,yy,40,'ok','filled',...
%     'MarkerEdgeColor',[0 0 0],...  %[0 0 0] [1 1 1]
%     'MarkerFaceColor',[1 1 1]);
s = scatter(xx,yy,40,'o','filled', 'MarkerEdgeColor',[0 0 0]);
colors = ones(length(xx), 3);
if exist('x0','var')
    map = [1 1 1; 0.9 0.9 0.9;
           1 0 0; 0.8 0   0;
           0 0 1; 0   0   0.8];
    for k=1:length(x0)
        colors(k,:) = map(x0(k),:);
    end
end


s.CData = colors;
title(tit)
hold off