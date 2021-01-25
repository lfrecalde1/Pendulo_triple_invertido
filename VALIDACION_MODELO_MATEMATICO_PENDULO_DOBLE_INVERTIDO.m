%% PROGRAMA PARA VALIDCION DEL MODELO MATEMATICO DE UN PENDULO DOBLE INVERTIDO
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
%% TIEMPOS DE SIMULACION 
ts=0.01;
t_final=20;
to=0;
t=[to:ts:t_final];
%% VALORES DEL SISTEMA PAPPER MARCELO
m_0=0.5;
m_1=0.4;
m_2=0.2;
l_1=0.2;
l_2=0.4;
g=9.8;
B_0=0.02;
B_1=0.02;
B_2=0.02;
I_1=(1/12)*m_1*l_1^2;
I_2=(1/12)*m_2*l_2^2;
N=(m_0+m_1+m_2)*g;
%% CONDICIONES INICIALES DEL SISTEMA
x(1)=0;
theta1(1)=0*pi/180;
theta2(1)=-1*pi/180;

x_p(1)=0;
theta1_p(1)=0*pi/180;
theta2_p(1)=0*pi/180;
%% ENTRADA DEL SISTEMA
%u=0.2*ones(1,length(t));

%% ESTADOS DESEADOS
xd=1*ones(1,length(t));
theta1d=0*ones(1,length(t));
theta2d=0*ones(1,length(t));

xd_p=0*ones(1,length(t));
theta1d_p=0*ones(1,length(t));
theta2d_p=0*ones(1,length(t));

xd_pp=0*ones(1,length(t));
theta1d_pp=0*ones(1,length(t));
theta2d_pp=0*ones(1,length(t));

K1=diag([1 1 1]);
K2=diag([1 1 1]);
K3=diag([5 5 5]);
K4=diag([1 1 1]);

%% BUCLE DE SIMULACION
for k=1:length(t)
    
   hd=[xd(k);theta1d(k);theta2d(k)];
   hdp=[xd_p(k);theta1d_p(k);theta2d_p(k)];
   hdpp=[xd_pp(k);theta1d_pp(k);theta2d_pp(k)];
   %% VECTORES DE VELOCIDADES
   hp=[x_p(k);theta1_p(k);theta2_p(k)];
   %% VECTORES DE POSICIONES
   h=[x(k);theta1(k);theta2(k)];
   %% ERRORES DEL SISTEMA
   
   he=hd-h;
   he_p=hdp-hp;
  
   %% VECTOR FRICCIONES
   f=[N;theta1_p(k);theta2_p(k)];
   M= [m_0 + m_1 + m_2, l_1*cos(theta1(k))*(m_1 + m_2), l_2*m_2*cos(theta2(k)); 
       l_1*cos(theta1(k))*(m_1 + m_2), I_1 + l_1^2*m_1 + l_1^2*m_2, l_1*l_2*m_2*cos(theta1(k) - theta2(k));
       l_2*m_2*cos(theta2(k)), l_1*l_2*m_2*cos(theta1(k) - theta2(k)), m_2*l_2^2 + I_2];
   
   C=[0, - l_1*m_1*theta1_p(k)*sin(theta1(k)) - l_1*m_2*theta1_p(k)*sin(theta1(k)), -l_2*m_2*theta2_p(k)*sin(theta2(k));
      0, 0, l_1*l_2*m_2*theta2_p(k)*sin(theta1(k) - theta2(k));
      0, -l_1*l_2*m_2*theta1_p(k)*sin(theta1(k) - theta2(k)),0];
  
   G=[0; -g*l_1*sin(theta1(k))*(m_1 + m_2); -g*l_2*m_2*sin(theta2(k))];
   
   E=[1;0;0];
   
   Fr=[-B_0*sign(x_p(k)), 0, 0; 0, -B_1, 0; 0, 0, -B_2];
   %% CONTROL DEL SISTEMA
   CONTROL=hdpp+K2*tanh(inv(K2)*K1*he)+K4*tanh(inv(K4)*K3*he_p);
%    u(k)=pinv(E)*(M*CONTROL+C*hp+G-Fr*f);
   u(k)=0;
   %% SITEMA SIMULACION
   hpp=pinv(M)*(E*u(k)+Fr*f-C*hp-G);
   %% INTEGRACION DE LOS ESTADOS DEL SISTEMA
   x_p(k+1)=x_p(k)+ts*hpp(1);
   theta1_p(k+1)=theta1_p(k)+ts*hpp(2);
   theta2_p(k+1)=theta2_p(k)+ts*hpp(3);
   
   x(k+1)=x(k)+ts*x_p(k);
   theta1(k+1)=theta1(k)+ts*theta1_p(k);
   theta2(k+1)=theta2(k)+ts*theta2_p(k);
   
   
end
figure
plot(t,x(1:length(t)),'r')
grid on
hold on
plot(t,theta1(1:length(t)),'b')
plot(t,theta2(1:length(t)),'g')
legend('posicion','angulo1','angulo2')