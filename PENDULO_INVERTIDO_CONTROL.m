%% VALIDACION DEL MODLEO MATEMATICO PENDULO SIMPLE INVERTIDO
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
%dbstop if error
%% TIEMPOS DE SIMULACION 
ts=0.01;
t_final=35;
to=0;
t=[to:ts:t_final];
%% VALORES DEL SISTEMA PAPPER MARCELO
m_0=0.48;
m_1=0.16;
l_1=0.5;
g=9.8;
B_0=3.83;
B_1=0.0218;
I_1=(1/12)*m_1*l_1^2;
%% CONDICIONES INICIALES DEL SISTEMA POSICIONES LINEALES Y ANGULARES
x(1)=2;
theta1(1)=165*pi/180;
%% CONDICIONES INICIALES DEL SISTEMA VELOCIDADES LINEALES Y ANGULARES
x_p(1)=0;
theta1_p(1)=0*pi/180;
%% ENTRADA DEL SISTEMA
u=0*ones(1,length(t));
%% BUCLE SE SOLUCION NUMERICA
I=diag([1 1]);
Z=diag([0 0]);
%% CONTROL 
A_LINEAL=[-4979/640, 327/4000, 0, -147/50;
  1149/80, -327/500, 0, 588/25;
  1, 0, 0, 0;
  0, 1, 0, 0];
B_LINEAL=[65/32; -15/4; 0; 0];
lambda=eig(A_LINEAL);
control=rank(ctrb(A_LINEAL,B_LINEAL));
Q=diag([1 1 10 10]);
R=1;
K=lqr(A_LINEAL,B_LINEAL,Q,R);
xr1=[-1*ones(1,500),-1*ones(1,500),-1*ones(1,500),-1*ones(1,500),0*ones(1,500),0*ones(1,500),0*ones(1,501)];
xr_p=0*ones(1,length(t));
phi=0*ones(1,length(t));
phi_p=0*ones(1,length(t));
wr = [xr_p; phi_p; xr1; phi]; % reference position
for k=1:length(t)
 h=[x_p(k);theta1_p(k);x(k);theta1(k)];
 he(:,k)=wr(:,k)-h;
 %% CALCULO DE LA ENERGIA DEL SISTEMA
 qp=[x_p(k);theta1_p(k)];
 M=[m_0 + m_1, l_1*m_1*cos(theta1(k)); 
     l_1*m_1*cos(theta1(k)), (13*l_1^2*m_1)/12];
 P=g*m_1*l_1*cos(theta1(k));
 Ep(k)=(1/2)*I_1*theta1_p(k)^2+m_1*g*l_1*(cos(theta1(k))-1);
 %% CONTROLADOR BASADOS EN FRECUENCIAS Y RETROALIMENTACION DE ESTADOS LINELAIZADOS
 error_angulo=phi(k)-Ang180(theta1(k));
 if(abs(error_angulo)>=20*pi/180)
    u(k)=9*sign(Ep(k)*theta1_p(k)*cos(theta1(k)));
 else
    u(k)=-K*(h - wr(:,k)); % control law 
 end
 X(:,k)=[x(k);x_p(k);theta1(k);theta1_p(k)];
 %% DINAMICA DEL SISTEMA EXPRESADO EN ESPACIOS DE ESTADOS
 A_1=[-(13*B_0)/(-12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1),(13*m_1*theta1_p(k)*sin(theta1(k))*l_1^2 + 12*B_1*cos(theta1(k)))/(l_1*(- 12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1));...
     (12*B_0*cos(theta1(k)))/(l_1*(-12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1)), -(6*(theta1_p(k)*sin(2*theta1(k))*l_1^2*m_1^2 + 2*B_1*m_1 + 2*B_1*m_0))/(l_1^2*m_1*(13*m_0 + 7*m_1 - 6*m_1*cos(2*theta1(k))))];
 A=[A_1,Z;I,Z];
 B_x=[13/(-12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1); 
     -(12*cos(theta1(k)))/(l_1*(- 12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1))];
 B=[B_x;0;0];
 G_1=[-(6*g*m_1*sin(2*theta1(k)))/(13*m_0 + 7*m_1 - 6*m_1*cos(2*theta1(k))); (12*g*sin(theta1(k))*(m_0 + m_1))/(l_1*(12*m_1*sin(theta1(k))^2 + 13*m_0 + m_1))];
 G=[G_1;0;0];
 
%% SIMULACION DEL SISTEMA 
 hp=A*h+B*u(k)+G;
 %% INTEGRACION NUMERICA PARA OBTENER LOS ESTADOS
 x_p(k+1)=x_p(k)+ts*hp(1);
 theta1_p(k+1)=theta1_p(k)+ts*hp(2);
 x(k+1)=x(k)+ts*hp(3);
 theta1(k+1)=theta1(k)+ts*hp(4);
 
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
legend('Angulo1')

figure(3)
for k=1:10:length(t)
    drawpend(X(:,k)',m_1,m_0,l_1);
    
end
figure
plot(t,he(1,1:length(t)),'r');
hold on;
grid on;
plot(t,he(2,1:length(t)),'b');
plot(t,he(3,1:length(t)),'g');
plot(t,he(4,1:length(t)),'k');
legend('error xp','error thp','error x','error th');
figure
plot(t,Ep,'--k')
grid on;
hold on
legend('Energia del sistema')