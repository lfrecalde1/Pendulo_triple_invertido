m=1.0;

m2=1.0;

L=1.0;
L2=1.0;

C=5;
g=9.8;
v0=0;
u=1;

theta0=pi/6;
alpha0=pi/2;

t0=0;
tf=100;

myopts = simset('MaxStep',0.02);
sim('pend', [t0 tf],myopts);


rball=.05;             % mass radius
rball1=.05;

x=L*sin(theta);
y=-L*cos(theta);

x1=x+L2*sin(alpha);
y1=y-L2*cos(alpha);

posx=x(1); 
posy=y(1);  
% Mass’s initial position
%Initialize figure, mass, and string

posx1=x1(1); 
posy1=y1(1);


fig=figure(2);
axs=axes('Parent',fig);

ball=rectangle('Position',[posx1-rball,posy1-rball,2*rball,2*rball],...
'Curvature',[1,1],...
'FaceColor','R',...
'Parent',axs);
rod=line([0 posx],[0 posy],'Marker','.','LineStyle','-');

ball1=rectangle('Position',[posx-rball1,posy-rball1,2*rball1,2*rball1],...
'Curvature',[1,1],...
'FaceColor','b',...
'Parent',axs);
rod1=line([posx posx1],[posy posy1],'Marker','.','LineStyle','-');



axis(axs,[-3,3,-3-rball,3]);
for j=2:length(time)
set(ball1,'Position',[x1(j)-rball1,y1(j)-rball1,2*rball1,2*rball1]);
set(rod1,'XData',[x(j) x1(j)],'YData',[y(j) y1(j)]);
set(ball,'Position',[x(j)-rball,y(j)-rball,2*rball,2*rball]);
set(rod,'XData',[0 x(j)],'YData',[0 y(j)]);
axis([-3,3,-3-rball,3])
pause(0.01);
end