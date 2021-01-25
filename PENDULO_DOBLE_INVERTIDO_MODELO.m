%% SIMULACION MODELO DINAMICA PENDULO INVERTIDO DOBLE
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
%dbstop if error
%% TIEMPOS DE SIMULACION 
ts=0.01;
t_final=35;
to=0;
t=[to:ts:t_final];
%% VALORES DEL SISTEMA PAPPER MARCELO
m_0=0.5;
m_1=0.1;
m_2=0.1;
l_1=0.25;
l_2=0.25;
g=9.8;
B_0=3.8;
B_1=0.0218;
B_2=0.0218;
I_1=(1/12)*m_1*l_1^2;
I_2=(1/12)*m_2*l_2^2;
%% CONDICIONES INICIALES DEL SISTEMA POSICIONES LINEALES Y ANGULARES
x(1)=0;
theta1(1)=0*pi/180;
theta2(1)=0*pi/180;
%% CONDICIONES INICIALES DEL SISTEMA VELOCIDADES LINEALES Y ANGULARES
x_p(1)=0;
theta1_p(1)=0*pi/180;
theta2_p(1)=0*pi/180;
%% ENTRADA DEL SISTEMA
u=0*ones(1,length(t));
%% CONTROL DEL SISTEMA
A_LINEAL=[-6878/919, 18312/114875, 1308/114875, 0, -16464/4595, -588/4595; 
          25536/919, -413328/114875,62784/22975,0,371616/4595, -28224/919; 
          1824/919, 62784/22975, -664464/114875, 0, -56448/919, 298704/4595; 
          1, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0; 
          0, 0, 1, 0, 0, 0];
B_LINEAL=[1810/919; -6720/919; -480/919; 0; 0; 0];
lambda=eig(A_LINEAL);
control=rank(ctrb(A_LINEAL,B_LINEAL));
Q=diag([1 1 1 10 10 10]);
R=1;
K=lqr(A_LINEAL,B_LINEAL,Q,R);
xr1=[0*ones(1,500),-1*ones(1,500),-1*ones(1,500),-1*ones(1,500),0*ones(1,500),0*ones(1,500),0*ones(1,501)];
xr_p=0*ones(1,length(t));
phi=0*ones(1,length(t));
phi_p=0*ones(1,length(t));
phi1=0*ones(1,length(t));
phi1_p=0*ones(1,length(t));
wr = [xr_p; phi_p;phi1_p; xr1; phi; phi1]; % reference position
for k=1:length(t)
    h=[x_p(k);theta1_p(k);theta2_p(k);x(k);theta1(k);theta2(k)];
    he(:,k)=wr(:,k)-h;
    u(k)=-K*(h - wr(:,k)); % control law
    if(t(k)>=30)
        u(k)=0;
    end
    %% ESTADOS PARA LA SIMULACION DEL SISTEMA 
    X(:,k)=[x(k);x_p(k);theta1(k);theta1_p(k);theta2(k);theta2_p(k)];
    %% VECTORES DE VELOCIDADES DLE SISTEMA
    hp=[x_p(k);theta1_p(k);theta2_p(k)];
    %% MATRICES DEL SISTEMA
    M=[m_0 + m_1 + m_2, l_1*cos(theta1(k))*(m_1 + m_2), l_2*m_2*cos(theta2(k));
       l_1*cos(theta1(k))*(m_1 + m_2), (l_1^2*(13*m_1 + 12*m_2))/12, l_1*l_2*m_2*cos(theta1(k) - theta2(k));
       l_2*m_2*cos(theta2(k)), l_1*l_2*m_2*cos(theta1(k) - theta2(k)), (13*l_2^2*m_2)/12];
   
    C=[B_0, - l_1*m_1*theta1_p(k)*sin(theta1(k)) - l_1*m_2*theta1_p(k)*sin(theta1(k)), -l_2*m_2*theta2_p(k)*sin(theta2(k));
       0, B_1, l_1*l_2*m_2*theta2_p(k)*sin(theta1(k) - theta2(k));
       0, -l_1*l_2*m_2*theta1_p(k)*sin(theta1(k) - theta2(k)), B_2];
   
   G=[0; -g*l_1*sin(theta1(k))*(m_1 + m_2); -g*l_2*m_2*sin(theta2(k))];
   E=[1;0;0];
%    u(k)=1.7*theta1_p(k);
   %% DINAMICA DEL SISTEMA SOLUCION NUMERICA
   hpp=inv(M)*(E*u(k)-C*hp-G);
   
   %% INTEGRACION NUMERICA PARA OBTENER LOS ESTADOS DEL SISTEMA
   x_p(k+1)=x_p(k)+ts*hpp(1);
   theta1_p(k+1)=theta1_p(k)+ts*hpp(2);
   theta2_p(k+1)=theta2_p(k)+ts*hpp(3);
   
   x(k+1)=x(k)+ts*x_p(k);
   theta1(k+1)=theta1(k)+ts*theta1_p(k);
   theta2(k+1)=theta2(k)+ts*theta2_p(k);
    
    
end
figure(3)
for k=1:10:length(t)  
    drawpend2(X(:,k)',m_1,m_2,m_0,l_1,l_2)
end
figure(1)
plot(t,x(1:length(t)),'r')
grid on
hold on
legend('posicion')
figure(2)
plot(t,theta1(1:length(t))*180/pi,'b')
grid on
hold on
plot(t,theta2(1:length(t))*180/pi,'g')
legend('Angulo1','Angulo2')
figure
plot(t,he(1,1:length(t)),'r');
hold on;
grid on;
plot(t,he(2,1:length(t)),'b');
plot(t,he(3,1:length(t)),'g');
plot(t,he(4,1:length(t)),'k');
plot(t,he(5,1:length(t)),'m');
plot(t,he(6,1:length(t)),'c');
legend('error xp','error th1p','error th2p','error x','error th1','error th2');
