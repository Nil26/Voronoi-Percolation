%clc;clear;
rng('shuffle');

% 32 42 48 64
for N = [32]
percolation_system_size(N);
end

% single_test();

function single_test()

x = gallery('uniformdata',[1 N],0);
y = gallery('uniformdata',[1 N],1);

X_init = [x;y];

x_p = [x-1,x-1,x-1,x,x,x,x+1,x+1,x+1];
y_p = [y-1,y,y+1,y-1,y,y+1,y-1,y,y+1];

X = [x_p;y_p];

p = 1;
vn = voronoi_neighbors(X',p);
vG = graph(full(vn));
[labels,~] = conncomp(vG);
end

function percolation_system_size(N)

generate_network_times = 10;

for net_time = 1:generate_network_times

fileID = fopen('delaunay_diff_network.txt','a');
fmt = '%5d %.4f\n';

x = gallery('uniformdata',[1 N],0);
y = gallery('uniformdata',[1 N],1);

X_init = [x;y];

x_p = [x-1,x-1,x-1,x,x,x,x+1,x+1,x+1];
y_p = [y-1,y,y+1,y-1,y,y+1,y-1,y,y+1];

X = [x_p;y_p];

average_times = 50;

for times = 1:average_times
p_first_time_reach_percolation(N,X_init,X,x,y,fileID,fmt);
end

end

fclose(fileID);
end

function p_first_time_reach_percolation(N,X_init,X,x,y,fileID,fmt)
% start from p0 to find the first time of p to reach percolation
span = 0;
p = 1/2;
delta = 1/2;
while(delta > 0.001 || span == 0)
delta = delta/2;
vn = voronoi_neighbors(X',p);
vG = graph(full(vn));
[labels,~] = conncomp(vG);
span = check_spanning_top_bottom(X_init,X,labels);
if span == 1
    p = p-delta;
else
    p = p+delta;
end

end

%fprintf(fileID,fmt,[N p]);

% draw_graph(vG,X,x,y,labels);
% title(sprintf("p = %.5f",p));

end

function draw_graph(vG,X,x,y,labels)

h = plot(vG,'LineWidth',2,'EdgeColor','blue');
h.MarkerSize = 10;
%h.NodeLabel = labels;
h.NodeCData = labels;
h.XData = X(1,:);
h.YData = X(2,:);

hold on; plot_voronoi(X,x,y);

end

function span = check_spanning_top_bottom(X_init,X,labels)

right_boundary = X_init(:,X_init(1,:) >= 0.9 & X_init(1,:) <= 1);
left_boundary = X_init(:,X_init(1,:) >= 0 & X_init(1,:) <= 0.1);

left_indices = logical(prod(ismember(X,left_boundary)));
right_indices = logical(prod(ismember(X,right_boundary)));

left_label = labels(left_indices);
right_label = labels(right_indices);

span = ~isempty(intersect(left_label,right_label));

% hold on;
% plot(right_boundary(1,:),right_boundary(2,:),'*');
% plot(left_boundary(1,:),left_boundary(2,:),'^');

end

function plot_voronoi(X,x,y)
[vx,vy] = voronoi(X(1,:),X(2,:));

bound_x = (min(x)+max(x))/2;
bound_y = (min(y)+max(y))/2;

plot(vx,vy,'black:','LineWidth',1);
axis equal;
xlim([-bound_x 3*bound_x]);
ylim([-bound_y 3*bound_y]);
hold on;
rectangle('Position',[0,0,1,1],...
         'LineWidth',2,'LineStyle','--','EdgeColor','r');

axis off;
set(gcf,'position',[500,500,500,500]);
set(gca,'LooseInset',get(gca,'TightInset'));
end


function vn = voronoi_neighbors(X,p)
% probability p to build the bond
  [n,~] = size(X);
%
%  V contains the Voronoi vertices, 
%  C contains indices of Voronoi vertices that form the (finite sides of the)
%  Voronoi cells.
%
  [~,C] = voronoin(X);
%
%  Two nodes are neighbors if they share an edge, that is, two Voronoi
%  vertices.
%
  vn = zeros(n,n);
  
  for i = 1 : n
    for j = i + 1 : n
      s = size(intersect(C{i},C{j}));
      if (1 < s(2)) && (rand < p)
        vn(i,j) = 1;
        vn(j,i) = 1;
      end
    end
  end
    
  return
end