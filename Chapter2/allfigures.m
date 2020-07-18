clear all;close all;clc;
figure(1);
axes1 = axes('Parent',figure(1));
set(axes1,'FontName','Times New Roman','FontSize',16);
hold(axes1,'on');



lambda=1 %/m^2 
diskRadius=100; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
diskArea=pi*diskRadius^2;

   randNumb_SBS_t1=poissrnd(lambda*diskArea);
   theta = rand(randNumb_SBS_t1,1)*(2*pi);
   r = diskRadius*sqrt(rand(randNumb_SBS_t1,1));
   x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y_1 =  r.*sin(theta);   %%%************************************************
   location_location = [x_1,y_1];
   
   
   plot(x_1,y_1,'.','markersize',15);
   
   xlim([-10,10]);
   ylim([-10,10]);
   axis('square');
      grid on;

   box on;
   saveas(figure(1),'PPP','eps');
   %%
   figure(2);
   axes1 = axes('Parent',figure(2));
   set(axes1,'FontName','Times New Roman','FontSize',16);
   hold(axes1,'on');
   m =50;
   sigma = 2;
   no_users_t1= poissrnd(m,1,1);
   user_pos=sigma*randn(sum(no_users_t1),2);
   plot(user_pos(:,1),user_pos(:,2),'.','markersize',15);
   %UE_location_all = [[0, 0];UE_location_al  xlim([-10,10]);
   ylim([-10,10]);
   xlim([-10,10]);
   axis('square');
      grid on;

   box on;
   saveas(figure(2),'PCP1','eps');
   %%
   figure(3);
   axes3 = axes('Parent',figure(3));
   set(axes3,'FontName','Times New Roman','FontSize',16);
   hold(axes3,'on');
   m =50;
   R = 5;
   no_users  = poissrnd(m,1,1);
   r = R*sqrt(rand(sum(no_users),1));
   theta = rand(sum(no_users),1)*(2*pi);
   x =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y =  r.*sin(theta);
   user_pos=[x,y];
   plot(user_pos(:,1),user_pos(:,2),'.','markersize',15);
   %UE_location_all = [[0, 0];UE_location_al  xlim([-10,10]);
   %%% Draw the circle
   ang=0:0.01:2*pi; 
   xp=R*cos(ang);
   yp=R*sin(ang);
   plot(xp,yp,'--','linewidth',2);
   
   
   grid on;
   ylim([-10,10]);
   xlim([-10,10]);
   
   axis('square');
   box on;
   
 annotation(figure(3),'textbox',...
    [0.552785714285713 0.450000000000001 0.0865 0.0619047619047619],...
    'String',{'$\mathtt{R}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off');

% Create arrow
annotation(figure(3),'arrow',[0.521428571428571 0.616071428571429],...
    [0.514285714285714 0.361904761904762],'LineWidth',2);
   
   
   saveas(figure(3),'PCP3','eps');
   
%% Example of PCP
figure(4);
axes4 = axes('Parent',figure(4));
set(axes4,'FontName','Times New Roman','FontSize',16);
hold(axes4,'on');



lambda=0.1 ;%/m^2 
diskRadius=20; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
diskArea=pi*diskRadius^2;

   randNumb_SBS_t1=poissrnd(lambda*diskArea);
   theta = rand(randNumb_SBS_t1,1)*(2*pi);
   r = diskRadius*sqrt(rand(randNumb_SBS_t1,1));
   x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y_1 =  r.*sin(theta);   %%%************************************************
   location_location = [x_1,y_1];
     
  l2 =  plot(x_1,y_1,'x','markersize',6,'linewidth',2);
   
   xlim([-10,10]);
   ylim([-10,10]);
   axis('square');
      grid on;

   for j=1:  randNumb_SBS_t1
       i = randi(3)
       switch i
           case 2
             l1 =  plot(x_1(j),y_1(j),'ko','LineWidth',3);
           case 3
               plot(-1+x_1(j),-1+y_1(j),'ko','LineWidth',3);
               plot(1+x_1(j),1+y_1(j),'ko','LineWidth',3);
               
              ang=0:0.01:2*pi; 
                xp= sqrt(2)*cos(ang);
                yp=sqrt(2)* sin(ang);
                 plot(x_1(j)+xp,y_1(j)+yp,'k:','linewidth',2);
               
       end
   end
   l = legend([l2,l1],'Parent point','Offspring point');
   set(l,'fontsize',16);   
   box on;
   saveas(figure(4),'PCP4','eps');

 %% MCP
 figure(5);
axes5 = axes('Parent',figure(5));
set(axes5,'FontName','Times New Roman','FontSize',16);
hold(axes5,'on');
lambda=0.1 ;%/m^2 

diskRadius=20; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
diskArea=pi*diskRadius^2;
randNumb_UE=poissrnd(lambda*diskArea);

theta = rand(randNumb_UE,1)*(2*pi);
r = diskRadius*sqrt(rand(randNumb_UE,1));
x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
y_1 =  r.*sin(theta);   %%%************************************************
UE_cc =[x_1, y_1];
UE_location_all=[];
total_user_count=0;

m=10;
no_users= poissrnd(m,randNumb_UE,1);
r = no_users; 
x = UE_cc;
t = r > 0;
a = cumsum(r(t));
b = zeros(1,a(end));
b(a - r(t) + 1) = 1;
x1 = UE_cc(t,:);
cc_location_rep = x1(cumsum(b),:);
R = 2;
theta = rand(sum(no_users),1)*(2*pi);
r = R*sqrt(rand(sum(no_users),1));
x =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
y =  r.*sin(theta);
user_pos=[x,y];
UE_location_all=cc_location_rep+ user_pos;

l2 =  plot(x_1,y_1,'x','markersize',6,'linewidth',2);
l1 = plot(UE_location_all(:,1),UE_location_all(:,2),'k.','markersize',10);

 ang=0:0.01:2*pi; 
 xp= R*cos(ang);
 yp=R* sin(ang);


for count  = 1: randNumb_UE
   plot(x_1(count)+xp,y_1(count)+yp,'k:');
end

  xlim([-10,10]);
   ylim([-10,10]);
   axis('square');
      grid on;
   l = legend([l2,l1],'Parent point','Offspring point');
   set(l,'Position',[0.211011911733519 0.109920638468531 0.307142850171243 0.129761901214009],...
    'FontSize',16);  
   box on;
   saveas(figure(5),'MCP','eps');  

 %% TCP
figure(6);
axes6 = axes('Parent',figure(6));
set(axes6,'FontName','Times New Roman','FontSize',16);
hold(axes6,'on');
lambda=0.1 ;%/m^2 

diskRadius=20; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
diskArea=pi*diskRadius^2;
randNumb_UE=poissrnd(lambda*diskArea);

theta = rand(randNumb_UE,1)*(2*pi);
r = diskRadius*sqrt(rand(randNumb_UE,1));
x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
y_1 =  r.*sin(theta);   %%%************************************************
UE_cc =[x_1, y_1];


m=10;
no_users= poissrnd(m,randNumb_UE,1);
r = no_users; 
x = UE_cc;
t = r > 0;
a = cumsum(r(t));
b = zeros(1,a(end));
b(a - r(t) + 1) = 1;
x1 = UE_cc(t,:);
cc_location_rep = x1(cumsum(b),:);
sigma = 0.8;
x =  sigma*(randn(sum(no_users),1));   %%%*****shifting origin to receiver location*******
y =  sigma*(randn(sum(no_users),1));
user_pos=[x,y];
UE_location_all=cc_location_rep+ user_pos;

l2 =  plot(x_1,y_1,'x','markersize',6,'linewidth',2);
l1 = plot(UE_location_all(:,1),UE_location_all(:,2),'k.','markersize',10);


  xlim([-10,10]);
   ylim([-10,10]);
   axis('square');
      grid on;
   l = legend([l2,l1],'Parent point','Offspring point');
   set(l,'Position',[0.211011911733519 0.109920638468531 0.307142850171243 0.129761901214009],...
    'FontSize',16);  
   box on;
   saveas(figure(6),'TCP','eps');  
 %% Voronoi PPP
 figure(7);
axes7 = axes('Parent',figure(7));
set(axes7,'FontName','Times New Roman','FontSize',16);
hold(axes7,'on');



lambda=0.1 %/m^2 
diskRadius=100; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
diskArea=pi*diskRadius^2;

   randNumb_SBS_t1=poissrnd(lambda*diskArea);
   theta = rand(randNumb_SBS_t1,1)*(2*pi);
   r = diskRadius*sqrt(rand(randNumb_SBS_t1,1));
   x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y_1 =  r.*sin(theta);   %%%************************************************
   location_location = [x_1,y_1];
   
   
   plot(x_1,y_1,'.','markersize',15);
   H= voronoi(x_1 ,y_1);
   H(2).LineWidth = 2;
  % H(2).Color = [1,1,1];
   xlim([-10,10]);
   ylim([-10,10]);
   axis('square');
      grid on;

   box on;
   saveas(figure(7),'PPP_voronoi','eps');
 %%
 
 
 