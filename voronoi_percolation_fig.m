clear;clc;

N = 16;

x = gallery('uniformdata',[1 N],0);
y = gallery('uniformdata',[1 N],1);

x_p = [x-1,x-1,x-1,x,x,x,x+1,x+1,x+1];
y_p = [y-1,y,y+1,y-1,y,y+1,y-1,y,y+1];

[vx,vy] = voronoi(x_p,y_p);

bound_x = (min(x)+max(x))/2;
bound_y = (min(y)+max(y))/2;

% figure;
% plot(vx,vy,'b-');
% % hold on;
% % plot(x,y,'r+');
% axis equal;

% Delaunay Triangulation
% DT = delaunayTriangulation(x',y');
% %figure;
% hold on;
% triplot(DT,'green');
% axis equal;

% xlim([-bound_x 3*bound_x]);
% ylim([-bound_y 3*bound_y]);
% hold on;
% rectangle('Position',[0,0,max(x)-min(x),max(x)-min(x)],...
%          'LineWidth',2,'LineStyle','--');

p1 = [vx(1,:);vy(1,:)];     % columns are "from" points
p2 = [vx(2,:);vy(2,:)];     % columns are "to" points

% use probability p
p = 0.2;

[vx_new,vy_new] = random_sample(p,p1,p2);

figure;
plot(vx_new,vy_new,'black-','LineWidth',2);
axis equal;
xlim([-bound_x 3*bound_x]);
ylim([-bound_y 3*bound_y]);
hold on;
rectangle('Position',[0,0,1,1],...
         'LineWidth',2,'LineStyle','--','EdgeColor','r');

axis off;
set(gcf,'position',[500,500,500,500]);
set(gca,'LooseInset',get(gca,'TightInset'));
% fig2 = gcf;
% F = getframe(gcf);
% BW = imbinarize(sum(F.cdata,3));

% % find the connected cluster
% [L,n] = bwlabel(BW);
% figure;
% mask = (L~=1);
% imagesc(L,'AlphaData',mask);
% axis equal;
% axis off;
% set(gcf,'position',[500,500,500,500]);
% set(gca,'LooseInset',get(gca,'TightInset'));
% hold on;
% rectangle('Position',[0,0,max(x)-min(x),max(x)-min(x)],...
%          'LineWidth',2,'LineStyle','--','EdgeColor','r');


% find someway to determine the giant connected cluster

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

function [vx_new,vy_new] = random_sample(p,p1,p2)
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