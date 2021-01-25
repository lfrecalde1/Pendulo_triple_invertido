%% SIMULACION MODELO DINAMICA PENDULO INVERTIDO DOBLE
%% BORAR VARIABLES DEL SISTEMA
clc,clear all,close all;
%% VARIABLES PARA LA COMUNICACION DEL SISTEMA
%% COMUNICACION MEMORIAS
loadlibrary('smClient64.dll','./smClient.h');
%% ABRIR MEMORIA COMPARTIDA
calllib('smClient64','openMemory','Sistema',2);
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
m_3=0.1;
l_1=0.25;
l_2=0.25;
l_3=0.25;
g=9.8;
B_0=3.8;
B_1=0.0218;
B_2=0.0218;
B_3=0.0218;
I_1=(1/12)*m_1*l_1^2;
I_2=(1/12)*m_2*l_2^2;
I_3=(1/12)*m_3*l_3^2;
%% CONDICIONES INICIALES DEL SISTEMA POSICIONES LINEALES Y ANGULARES
x(1)=0;
theta1(1)=0*pi/180;
theta2(1)=0*pi/180;
theta3(1)=2*pi/180;
%% CONDICIONES INICIALES DEL SISTEMA VELOCIDADES LINEALES Y ANGULARES
x_p(1)=0;
theta1_p(1)=0*pi/180;
theta2_p(1)=0*pi/180;
theta3_p(1)=0*pi/180;
%% CONTROL DEL SISTEMA
A_LINEAL=[[ -38/5,   109/625,        0,        0, 0, -147/25,      0,      0]
[ 152/5, -2616/625,  436/125,        0, 0, 3528/25, -392/5,      0]
[     0,   436/125, -872/125,  436/125, 0,  -588/5,  784/5, -196/5]
[     0,         0,  436/125, -872/125, 0,       0, -392/5,  392/5]
[     1,         0,        0,        0, 0,       0,      0,      0]
[     0,         1,        0,        0, 0,       0,      0,      0]
[     0,         0,        1,        0, 0,       0,      0,      0]
[     0,         0,        0,        1, 0,       0,      0,      0]];

B_LINEAL=[2;
  -8;
   0;
   0;
   0;
   0;
   0;
   0];
u=0*ones(1,length(t));
lambda=eig(A_LINEAL);
control=rank(ctrb(A_LINEAL,B_LINEAL));
Q=diag([1 1 1 1 10 10 10 10]);
R=1;
K=lqr(A_LINEAL,B_LINEAL,Q,R);
xr1=[0*ones(1,500),-1*ones(1,500),-1*ones(1,500),0*ones(1,500),0*ones(1,500),1*ones(1,500),1*ones(1,501)];
xr_p=0*ones(1,length(t));
phi=0*ones(1,length(t));
phi_p=0*ones(1,length(t));
phi1=0*ones(1,length(t));
phi1_p=0*ones(1,length(t));
phi2=0*ones(1,length(t));
phi2_p=0*ones(1,length(t));
wr = [xr_p; phi_p;phi1_p ;phi2_p; xr1; phi; phi1;phi2]; % posicoon referencia

posicion(1)=x(1);
angulo_1(1)=Ang180(theta1(1));
angulo_2(1)=Ang180(theta2(1));
angulo_3(1)=Ang180(theta3(1));
%% GRAFICAS EN TIEMPO REAL
% figure(1)
% axis equal;
% h1=plot(t(1),posicion(1),'--r'); 
% hold on;
% h2=plot(t(1),angulo_1(1),'--g');
% grid on;
% h3=plot(t(1),angulo_2(1),'--b');
% grid on;
% h4=plot(t(1),angulo_3(1),'--c');
% hold on;
% title('GRAFICAS SISTEMA');
% xlabel('tiempo [s]'); ylabel('estados');
for k=1:length(t)
    tic;
    %% ACTUALIZACION PARA LAS GRADICAS EN TIEMPO REAL
%     drawnow;
%     delete(h1);
%     delete(h2);
%     delete(h3);
%     delete(h4);
    %% CONVERSION PARA LA COMUNICACION CON UNITY
    posicion(k)=x(k);
    angulo_1(k)=Ang180(theta1(k));
    angulo_2(k)=Ang180(theta2(k));
    angulo_3(k)=Ang180(theta3(k));
    %% ESCRITURA DE DE LOS ESTADOS DEL SISTEMA
    calllib('smClient64','setFloat','Sistema',0,posicion(k));
    calllib('smClient64','setFloat','Sistema',1,angulo_1(k));
    calllib('smClient64','setFloat','Sistema',2,angulo_2(k));
    calllib('smClient64','setFloat','Sistema',3,angulo_3(k));
    
    h=[x_p(k);theta1_p(k);theta2_p(k);theta3_p(k);x(k);theta1(k);theta2(k);theta3(k)];
    he(:,k)=wr(:,k)-h;
    u(k)=-K*(h - wr(:,k)); % control law
    if(t(k)>=25)
        u(k)=0;
    end
    %% ESTADOS PARA LA SIMULACION DEL SISTEMA 
    X(:,k)=[x(k);x_p(k);theta1(k);theta1_p(k);theta2(k);theta2_p(k);theta3(k);theta3_p(k)];
    %% VECTORES DE VELOCIDADES DLE SISTEMA
    hp=[x_p(k);theta1_p(k);theta2_p(k);theta3_p(k)];
    %% MATRICES DEL SISTEMA
    M=[m_0 + m_1 + m_2 + m_3, l_1*cos(theta1(k))*(m_1 + m_2 + m_3), l_2*cos(theta2(k))*(m_2 + m_3), l_3*m_3*cos(theta3(k));
       l_1*cos(theta1(k))*(m_1 + m_2 + m_3), l_1^2*(m_1 + m_2 + m_3), l_1*l_2*cos(theta1(k) - theta2(k))*(m_2 + m_3), l_1*l_3*m_3*cos(theta1(k) - theta3(k));
       l_2*cos(theta2(k))*(m_2 + m_3), l_1*l_2*cos(theta1(k) - theta2(k))*(m_2 + m_3), l_2^2*(m_2 + m_3), l_2*l_3*m_3*cos(theta2(k) - theta3(k));
       l_3*m_3*cos(theta3(k)), l_1*l_3*m_3*cos(theta1(k) - theta3(k)), l_2*l_3*m_3*cos(theta2(k) - theta3(k)), l_3^2*m_3];
   
   C=[B_0, - l_1*m_1*theta1_p(k)*sin(theta1(k)) - l_1*m_2*theta1_p(k)*sin(theta1(k)) - l_1*m_3*theta1_p(k)*sin(theta1(k)), - l_2*m_2*theta2_p(k)*sin(theta2(k)) - l_2*m_3*theta2_p(k)*sin(theta2(k)), -l_3*m_3*theta3_p(k)*sin(theta3(k));
      0, B_1, l_1*l_2*m_2*theta2_p(k)*sin(theta1(k) - theta2(k)) + l_1*l_2*m_3*theta2_p(k)*sin(theta1(k) - theta2(k)), l_1*l_3*m_3*theta3_p(k)*sin(theta1(k) - theta3(k));
      0, - l_1*l_2*m_2*theta1_p(k)*sin(theta1(k) - theta2(k)) - l_1*l_2*m_3*theta1_p(k)*sin(theta1(k) - theta2(k)), B_2, l_2*l_3*m_3*theta3_p(k)*sin(theta2(k) - theta3(k));
      0, -l_1*l_3*m_3*theta1_p(k)*sin(theta1(k) - theta3(k)), -l_2*l_3*m_3*theta2_p(k)*sin(theta2(k) - theta3(k)), B_3];
   
   G=[0;
      -g*l_1*sin(theta1(k))*(m_1 + m_2 + m_3);
      -g*l_2*sin(theta2(k))*(m_2 + m_3);
      -g*l_3*m_3*sin(theta3(k))];
  
   E=[1;0;0;0];

   %% DINAMICA DEL SISTEMA SOLUCION NUMERICA
   hpp=inv(M)*(E*u(k)-C*hp-G);
   
   %% INTEGRACION NUMERICA PARA OBTENER LOS ESTADOS DEL SISTEMA
   x_p(k+1)=x_p(k)+ts*hpp(1);
   theta1_p(k+1)=theta1_p(k)+ts*hpp(2);
   theta2_p(k+1)=theta2_p(k)+ts*hpp(3);
   theta3_p(k+1)=theta3_p(k)+ts*hpp(4);
   
   x(k+1)=x(k)+ts*x_p(k);
   theta1(k+1)=theta1(k)+ts*theta1_p(k);
   theta2(k+1)=theta2(k)+ts*theta2_p(k);
   theta3(k+1)=theta3(k)+ts*theta3_p(k);
   %% GRAFICASD E ESTADOS DEL SISTEMA
%    h1=plot(t(1:k),posicion(1:k),"--r");
%    hold on;
%    h2=plot(t(1:k),angulo_1(1:k),"--g");
%    hold on;
%    h3=plot(t(1:k),angulo_2(1:k),"--b");
%    hold on;
%    h4=plot(t(1:k),angulo_3(1:k),"--c");
%    hold on; 
%    legend('Posicion','Angulo 1','Angulo 2','Angulo 3');
%    grid on;
   while(toc<ts)
   end
   toc
end
%% LIBERAR MEMORIA COMPARTIDA
calllib('smClient64','freeViews')
unloadlibrary smClient64
figure(2)
for k=1:10:length(t)  
    drawpend2(X(:,k)',m_1,m_2,m_3,m_0,l_1,l_2,l_3)
end
figure(3)
plot(t,x(1:length(t)),'r')
grid on
hold on
legend('posicion')
figure(4)
plot(t,theta1(1:length(t))*180/pi,'b')
grid on
hold on
plot(t,theta2(1:length(t))*180/pi,'g')
plot(t,theta3(1:length(t))*180/pi,'m')
legend('Angulo1','Angulo2','Angulo3')
figure
plot(t,he(1,1:length(t)),'r');
hold on;
grid on;
plot(t,he(2,1:length(t)),'b');
plot(t,he(3,1:length(t)),'g');
plot(t,he(4,1:length(t)),'k');
plot(t,he(5,1:length(t)),'m');
plot(t,he(6,1:length(t)),'c');
plot(t,he(7,1:length(t)),'m');
plot(t,he(8,1:length(t)),'c');
legend('error xp','error th1p','error th2p','error th3p','error x','error th1','error th2','error th3');
