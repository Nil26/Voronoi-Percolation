clear;clc;
rng('shuffle');

% for voronoi

% 32 42 48 64
for N = [256]
percolation_system_size(N);
end

function percolation_system_size(N)

generate_network_times = 10;

for net_time = 1:generate_network_times

fileID = fopen('voronoi_diff_network.txt','a');
fmt = '%5d %.4f\n';

x = gallery('uniformdata',[1 N],0);
y = gallery('uniformdata',[1 N],1);

average_times = 50;

for times = 1:average_times
p_first_time_reach_percolation(x,y,N,fileID,fmt);
end

end

fclose(fileID);
end

function p_first_time_reach_percolation(x,y,N,fileID,fmt)
% start from p0 to find the first time of p to reach percolation
span = 0;
p = 1/2;
delta = 1/2;
while(delta > 0.001 || span == 0)
delta = delta/2;

x_p = [x-1,x-1,x-1,x,x,x,x+1,x+1,x+1];
y_p = [y-1,y,y+1,y-1,y,y+1,y-1,y,y+1];

[vx,vy] = voronoi(x_p,y_p);

p1 = [vx(1,:);vy(1,:)];     % columns are "from" points
p2 = [vx(2,:);vy(2,:)];     % columns are "to" points

%[vx_new,vy_new] = site_percolation(p,p1,p2);
% for bond percolation
vx_new = vx;
vy_new = vy;

p1_new = [vx_new(1,:);vy_new(1,:)];     % columns are "from" points
p2_new = [vx_new(2,:);vy_new(2,:)];     % columns are "to" points

all_vertices = union(p1_new',p2_new','rows')';
% hold on;
% plot(all_vertices(1,:),all_vertices(2,:),'^','MarkerSize',10);

% the adjacency matrix
n_vertices = size(all_vertices,2);
vn = zeros(n_vertices,n_vertices);

n_bonds = size(vx_new,2);
for e = 1:n_bonds
    pi = [vx_new(1,e);vy_new(1,e)];
    pj = [vx_new(2,e);vy_new(2,e)];
    [~,i_index] = ismember(pi',all_vertices','rows');
    [~,j_index] = ismember(pj',all_vertices','rows');
    if(rand < p)
    vn(i_index,j_index) = 1;
    vn(j_index,i_index) = 1;
    end
end

% build the graph and find components
vG=graph(vn);
[labels,~] = conncomp(vG);

% check spanning
span = check_spanning_top_bottom(all_vertices,labels);

if span == 1
    p = p-delta;
else
    p = p+delta;
end

end

fprintf(fileID,fmt,[N p]);

% plotting_all(x,y,vx_new,vy_new,vG,labels,all_vertices);
% title(sprintf("p = %.5f",p));

end


function plotting_all(x,y,vx_new,vy_new,vG,labels,all_vertices)

bound_x = (min(x)+max(x))/2;
bound_y = (min(y)+max(y))/2;

figure;
plot(vx_new,vy_new,'black:','LineWidth',1);
axis equal;
xlim([-bound_x 3*bound_x]);
ylim([-bound_y 3*bound_y]);
hold on;
rectangle('Position',[0,0,1,1],...
         'LineWidth',2,'LineStyle','--','EdgeColor','r');

axis off;
set(gcf,'position',[500,500,500,500]);
set(gca,'LooseInset',get(gca,'TightInset'));

hold on;
h = plot(vG,'LineWidth',2,'EdgeColor','black');
h.MarkerSize = 10;
%h.NodeLabel = labels;
h.NodeCData = labels;
h.XData = all_vertices(1,:);
h.YData = all_vertices(2,:);

end

function span = check_spanning_top_bottom(X,labels)

right_boundary = X(:,X(1,:) >= 0.9 & X(1,:) <= 1 & X(2,:) <= 1 & X(2,:) >= 0);
left_boundary = X(:,X(1,:) >= 0 & X(1,:) <= 0.1 & X(2,:) <= 1 & X(2,:) >= 0);

left_indices = logical(prod(ismember(X,left_boundary)));
right_indices = logical(prod(ismember(X,right_boundary)));

left_label = labels(left_indices);
right_label = labels(right_indices);

span = ~isempty(intersect(left_label,right_label));

% hold on;
% plot(right_boundary(1,:),right_boundary(2,:),'*');
% plot(left_boundary(1,:),left_boundary(2,:),'^');

end

function [p1_new,p2_new] = pickout_node(node,p1,p2)
%[~,ind1]=ismembertol(node,p1','ByRows',true);
% pick out the node from the diagram
ind1 = find(all(bsxfun(@eq, p1', node),2)); %// only where all line matches
p1(:,ind1) = [];
p2(:,ind1) = [];
ind2 = find(all(bsxfun(@eq, p2', node),2)); %// only where all line matches
p1(:,ind2) = [];
p2(:,ind2) = [];

p1_new = p1;
p2_new = p2;
end

function [vx_new,vy_new] = site_percolation(p,p1,p2)
num_node = length(p1);
rand_num = rand(1,num_node);
index_found = find(rand_num > p); % index_found means the nodes which are un-occupied
p1_original = p1;
for n = index_found
node = p1_original(:,n)';
[p1,p2] = pickout_node(node,p1,p2);
end

vx_new = [p1(1,:);p2(1,:)];
vy_new = [p1(2,:);p2(2,:)];
end