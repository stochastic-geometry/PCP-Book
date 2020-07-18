clear all; close all; clc;clearvars;

L=1
bounds = [0, L, L, 0; 0, 0, L, L];  % Bounding region

points1=[0.5    0.2    0.25   0.7; 0.6    0.8    0.3   0.2];

points1PP=[0.3    0.1    0.15   0.9, .6; 0.6    0.8    0.3   0.2,0.7];

sigma=1;

B_x2=[];
B_y2=[];
N_c2=size(points1PP,2)
for count1=1:N_c2
  no_ue=poissrnd(20)
  for i=1:no_ue
%       Rx=RR.*sqrt(rand(1,1));           %%%matern
%       thetaX=rand(1,1).*2*pi;           %%%matern
%       Rx=[Rx*cos(thetaX) Rx*sin(thetaX)];  
      B_x2=[B_x2 points1PP(1, count1)+sigma*randn];
      B_y2=[B_y2 points1PP(2, count1)+sigma*randn];
  end
end

N_s3=size(B_x2,2);
N_m1 = 4;                           % No of macrocells
N3 = N_m1 + N_s3; 

poinst2= [B_x2;B_y2]

points_MAS1 =[points1,poinst2];            % BS locations
W_m = 5.6;                            % Weight of macro
W_s = 1;                            % Weight of small cell

points_MAS1(3, 1:N_m1 ) = W_m;                % Weight of a Macrocells
points_MAS1(3, N_m1+1:N3) = W_s;              % Weight of a Small cells
%%  Model 3
U_xx2=[];
U_yy2=[];
sigma_u = 0.03;
for count1=1:size(poinst2,2)
  no_ue=poissrnd(20)
  for i=1:no_ue
%       Rx=RR.*sqrt(rand(1,1));           %%%matern
%       thetaX=rand(1,1).*2*pi;           %%%matern
%       Rx=[Rx*cos(thetaX) Rx*sin(thetaX)];  
      U_xx2=[U_xx2 poinst2(1, count1)+sigma_u*randn];
      U_yy2=[U_yy2 poinst2(2, count1)+sigma_u*randn];
  end
end

a_color=[1,0,0]

figure(1)
hold on
regionData = mwvoronoi(bounds, points_MAS1);

% Draw the required regions. The code of drawRegions is customized.
drawRegions(bounds, regionData, 'linewidth', 2, 'Color', 'b', 'LineStyle', '-');



plot(points_MAS1(1, N_m1+1:N3), points_MAS1(2,N_m1+1:N3), 'ko','MarkerFaceColor',[0,0,0]);
plot(U_xx2,U_yy2,'o','MarkerFaceColor',a_color,'MarkerEdgeColor', a_color,'Markersize',3)

aa=.8;
plot(points_MAS1(1,1:N_m1), points_MAS1(2,1:N_m1),'s','MarkerFaceColor',[aa,aa,aa],'MarkerEdgeColor',[0 0 0], 'Markersize',15,'Linewidth',2);
plot(points_MAS1(1,1:N_m1), points_MAS1(2,1:N_m1), 'ks','MarkerFaceColor',[0,0,0]);

%scatter(points1PP(1,:),points1PP(2,:),'ko','filed', 'SizeData',40)
xlim([min(bounds(1,:)), max(bounds(1,:))]);
ylim([min(bounds(2,:)), max(bounds(2,:))]);
set(gca, 'XTick', [], 'YTick', []);     % If Ticks on the axis are not desired
axis square
saveas(figure(1), 'Model2_max_power','eps') 

%% Model 4
U_PP=L * rand(2, 40*5);

figure(2)

U_xx2=[];
U_yy2=[];
sigma_u = 1;
for count1=1:size(poinst2,2)
  no_ue=poissrnd(20)
  for i=1:no_ue
%       Rx=RR.*sqrt(rand(1,1));           %%%matern
%       thetaX=rand(1,1).*2*pi;           %%%matern
%       Rx=[Rx*cos(thetaX) Rx*sin(thetaX)];  
      U_xx2=[U_xx2 poinst2(1, count1)+sigma_u*randn];
      U_yy2=[U_yy2 poinst2(2, count1)+sigma_u*randn];
  end
end



hold on
regionData = mwvoronoi(bounds, points_MAS1);

% Draw the required regions. The code of drawRegions is customized.
drawRegions(bounds, regionData, 'linewidth', 2, 'Color', 'b', 'LineStyle', '-');


plot(points_MAS1(1, N_m1+1:N3), points_MAS1(2,N_m1+1:N3), 'ko','MarkerFaceColor',[0,0,0]);
plot(U_xx2,U_yy2,'o','MarkerFaceColor',a_color,'MarkerEdgeColor', a_color,'Markersize',3)

aa=.8;
plot(points_MAS1(1,1:N_m1), points_MAS1(2,1:N_m1),'s','MarkerFaceColor',[aa,aa,aa],'MarkerEdgeColor',[0 0 0], 'Markersize',15,'Linewidth',2);
plot(points_MAS1(1,1:N_m1), points_MAS1(2,1:N_m1), 'ks','MarkerFaceColor',[0,0,0]);

%scatter(points1PP(1,:),points1PP(2,:),'ko','filed', 'SizeData',40)
xlim([min(bounds(1,:)), max(bounds(1,:))]);
ylim([min(bounds(2,:)), max(bounds(2,:))]);
set(gca, 'XTick', [], 'YTick', []);     % If Ticks on the axis are not desired
axis square
saveas(figure(2), 'Model1_max_power','eps') 

