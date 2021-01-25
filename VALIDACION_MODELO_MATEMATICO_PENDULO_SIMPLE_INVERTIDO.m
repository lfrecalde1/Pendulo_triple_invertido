%% VALIDACION DEL MODLEO MATEMATICO PENDULO SIMPLE INVERTIDO
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
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
for k=1:length(t)
   X(:,k)=[x(k);x_p(k);theta1(k);theta1_p(k)];
   h=[x(k);theta1(k)];
   hp=[x_p(k);theta1_p(k)];
   M=[m_0 + m_1, l_1*m_1*cos(theta1(k)); l_1*m_1*cos(theta1(k)), m_1*l_1^2 + I_1];
   C=[0, -l_1*m_1*theta1_p(k)*sin(theta1(k)); 0, 0];
   G=[0; -g*l_1*m_1*sin(theta1(k))];
   E=[1;0];
   Fr=[-B_0, 0; 0, -B_1];
   u(k)=2.542*theta1_p(k);
   hpp=inv(M)*(E*u(k)+Fr*hp-C*hp-G);
   
   x_p(k+1)=x_p(k)+ts*hpp(1);
   theta1_p(k+1)=theta1_p(k)+ts*hpp(2);
   
   x(k+1)=x(k)+ts*x_p(k);
   theta1(k+1)=theta1(k)+ts*theta1_p(k);
   
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