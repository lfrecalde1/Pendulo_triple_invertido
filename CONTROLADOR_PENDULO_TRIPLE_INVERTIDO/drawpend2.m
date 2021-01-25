function drawpend2(state,m_1,m_2,m_3,m_0,l_1,l_2,l_3)
x = state(1);
th = state(3);
th1=state(5);
th2=state(7);
% dimensions
W = 1*sqrt(m_0/5);  % cart width
H = .5*sqrt(m_0/5); % cart height
wr = .1;          % wheel radius
mr =0.3*sqrt(m_1);  % mass radius
mr1=0.3*sqrt(m_2);
mr2=0.3*sqrt(m_3);
% positions
y = wr/2+H/2; % cart vertical position
pendx = x + l_1*sin(th);
pendy = y + l_1*cos(th);

pendx1 = x + l_1*sin(th)+l_2*sin(th1);
pendy1 = y + l_1*cos(th)+l_2*cos(th1);

pendx2 = x + l_1*sin(th)+l_2*sin(th1)+l_3*sin(th2);
pendy2 = y + l_1*cos(th)+l_2*cos(th1)+l_3*cos(th2);


plot([-10 10],[0 0],'k','LineWidth',2), hold on

rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[.5 0.5 1],'LineWidth',1.5); % Draw cart
rectangle('Position',[x-.9*W/2,0,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5); % Draw wheel
rectangle('Position',[x+.9*W/2-wr,0,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5); % Draw wheel

plot([x pendx],[y pendy],'k','LineWidth',2); % Draw pendulum
plot([pendx pendx1],[pendy pendy1],'k','LineWidth',2); % Draw pendulum
plot([pendx1 pendx2],[pendy1 pendy2],'k','LineWidth',2); % Draw pendulum

rectangle('Position',[pendx-mr/2,pendy-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5);
rectangle('Position',[pendx1-mr1/2,pendy1-mr1/2,mr1,mr1],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5);
rectangle('Position',[pendx2-mr2/2,pendy2-mr2/2,mr2,mr2],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5);

axis([-5 5 -1 1]);, axis equal
grid on;
set(gcf,'Position',[100 100 1000 800])
drawnow, hold off