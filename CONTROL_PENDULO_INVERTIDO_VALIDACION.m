%% VALIDACION DEL MODLEO MATEMATICO PENDULO SIMPLE INVERTIDO
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
dbstop if error
%% TIEMPOS DE SIMULACION 
ts=0.01;
t_final=5;
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
x(1)=0;
theta1(1)=165*pi/180;
%% CONDICIONES INICIALES DEL SISTEMA VELOCIDADES LINEALES Y ANGULARES
x_p(1)=0;
theta1_p(1)=0*pi/180;
%% ENTRADA DEL SISTEMA
u=0*ones(1,length(t));
%% BUCLE SE SOLUCION NUMERICA
I=diag([1 1]);
Z=diag([0 0]);

for k=1:length(t)
 h=[x_p(k);theta1_p(k);x(k);theta1(k)];
 
 X(:,k)=[x(k);x_p(k);theta1(k);theta1_p(k)];
 u(k)=2.541*theta1_p(k);
 
hp=[(13*u(k)*l_1 - 13*B_0*l_1*x_p(k) + 12*B_1*theta1_p(k)*cos(theta1(k)) + 13*l_1^2*m_1*theta1_p(k)^2*sin(theta1(k)) - 12*g*l_1*m_1*cos(theta1(k))*sin(theta1(k)))/(l_1*(- 12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1));
    -(12*B_1*m_0*theta1_p(k) + 12*B_1*m_1*theta1_p(k) + 12*u(k)*l_1*m_1*cos(theta1(k)) - 12*g*l_1*m_1^2*sin(theta1(k)) + 12*l_1^2*m_1^2*theta1_p(k)^2*cos(theta1(k))*sin(theta1(k)) - 12*B_0*l_1*m_1*x_p(k)*cos(theta1(k)) - 12*g*l_1*m_0*m_1*sin(theta1(k)))/(l_1^2*m_1*(- 12*m_1*cos(theta1(k))^2 + 13*m_0 + 13*m_1));
     x_p(k); 
     theta1_p(k)]
 
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
for k=1:5:length(t)
    drawpend(X(:,k)',m_1,m_0,l_1);
end