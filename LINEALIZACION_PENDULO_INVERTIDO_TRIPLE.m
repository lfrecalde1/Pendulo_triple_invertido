clc;
clear all;
close all;
syms F g I_1 I_2 l_1 l_2 l_3 m_0 m_1 m_2 m_3 x_o x_op x_opp theta1 theta_1p theta_1pp theta2 theta_2p theta_2pp theta3 theta_3p theta_3pp B_0 B_1 B_2 B_3 as real
%Calculo de Posciones y Cinetica
hx1=x_o;
hy1=0;

h1=[hx1;hy1]  ;

hx2=x_o+l_1*sin(theta1);
hy2=l_1*cos(theta1);

h2=[hx2;hy2];

hx3=x_o+l_1*sin(theta1)+l_2*sin(theta2);
hy3=l_1*cos(theta1)+l_2*cos(theta2);

h3=[hx3;hy3];

hx4=x_o+l_1*sin(theta1)+l_2*sin(theta2)+l_3*sin(theta3);
hy4=l_1*cos(theta1)+l_2*cos(theta2)+l_3*cos(theta3);

h4=[hx4;hy4];


hx1p=simplify(diff(hx1,x_o)*x_op+diff(hx1,theta1)*theta_1p+diff(hx1,theta2)*theta_2p+diff(hx1,theta3)*theta_3p);
hy1p=simplify(diff(hy1,x_o)*x_op+diff(hy1,theta1)*theta_1p+diff(hy1,theta2)*theta_2p+diff(hy1,theta3)*theta_3p);
h1p=[hx1p;hy1p];
% 
hx2p=simplify(diff(hx2,x_o)*x_op+diff(hx2,theta1)*theta_1p+diff(hx2,theta2)*theta_2p+diff(hx2,theta3)*theta_3p);
hy2p=simplify(diff(hy2,x_o)*x_op+diff(hy2,theta1)*theta_1p+diff(hy2,theta2)*theta_2p+diff(hy2,theta3)*theta_3p);
h2p=[hx2p;hy2p];
% 
% 
hx3p=simplify(diff(hx3,x_o)*x_op+diff(hx3,theta1)*theta_1p+diff(hx3,theta2)*theta_2p+diff(hx3,theta3)*theta_3p);
hy3p=simplify(diff(hy3,x_o)*x_op+diff(hy3,theta1)*theta_1p+diff(hy3,theta2)*theta_2p+diff(hy3,theta3)*theta_3p);
h3p=[hx3p;hy3p];
%
hx4p=simplify(diff(hx4,x_o)*x_op+diff(hx4,theta1)*theta_1p+diff(hx4,theta2)*theta_2p+diff(hx4,theta3)*theta_3p);
hy4p=simplify(diff(hy4,x_o)*x_op+diff(hy4,theta1)*theta_1p+diff(hy4,theta2)*theta_2p+diff(hy4,theta3)*theta_3p);
h4p=[hx4p;hy4p];
%Calculo de Energia Cinetica Total
K=simplify(1/2*h1p'*[m_0,0;0,m_0]*h1p+1/2*h2p'*[m_1,0;0,m_1]*h2p+1/2*h3p'*[m_2,0;0,m_2]*h3p+1/2*h4p'*[m_3,0;0,m_3]*h4p)
P=simplify(m_1*g*hy2+m_2*g*hy3+m_3*g*hy4)
L=simplify(K-P);
% 
%Calculo del Lagrangiano 
lxp=diff(L,x_op);
T1=diff(lxp,x_o)*x_op+diff(lxp,theta1)*theta_1p+diff(lxp,theta2)*theta_2p+diff(lxp,theta3)*theta_3p+diff(lxp,x_op)*x_opp+diff(lxp,theta_1p)*theta_1pp+diff(lxp,theta_2p)*theta_2pp+diff(lxp,theta_3p)*theta_3pp-diff(L,x_o);
T1=simplify(T1);
% 
ltheta1p=diff(L,theta_1p);
T2=simplify(diff(ltheta1p,x_o)*x_op+diff(ltheta1p,theta1)*theta_1p+diff(ltheta1p,theta2)*theta_2p+diff(ltheta1p,theta3)*theta_3p+diff(ltheta1p,x_op)*x_opp+diff(ltheta1p,theta_1p)*theta_1pp+diff(ltheta1p,theta_2p)*theta_2pp+diff(ltheta1p,theta_3p)*theta_3pp-diff(L,theta1));
T2=simplify(T2);
% 
ltheta2p=diff(L,theta_2p);
T3=simplify(diff(ltheta2p,x_o)*x_op+diff(ltheta2p,theta1)*theta_1p+diff(ltheta2p,theta2)*theta_2p+diff(ltheta2p,theta3)*theta_3p+diff(ltheta2p,x_op)*x_opp+diff(ltheta2p,theta_1p)*theta_1pp+diff(ltheta2p,theta_2p)*theta_2pp+diff(ltheta2p,theta_3p)*theta_3pp-diff(L,theta2));
T3=simplify(T3);
%
ltheta3p=diff(L,theta_3p);
T4=simplify(diff(ltheta3p,x_o)*x_op+diff(ltheta3p,theta1)*theta_1p+diff(ltheta3p,theta2)*theta_2p+diff(ltheta3p,theta3)*theta_3p+diff(ltheta3p,x_op)*x_opp+diff(ltheta3p,theta_1p)*theta_1pp+diff(ltheta3p,theta_2p)*theta_2pp+diff(ltheta3p,theta_3p)*theta_3pp-diff(L,theta3));
T4=simplify(T4);

%Calculo de matriz M
M=simplify(jacobian([T1,T2,T3,T4],[x_opp, theta_1pp,theta_2pp,theta_3pp]));
Mqpp=M*[x_opp;theta_1pp;theta_2pp;theta_3pp];
%Calculo de la matriz C y G
T=[T1;T2;T3;T4];
CG=simplify(T-Mqpp);
CG_0(g)=CG;
Cqp=simplify(CG_0(0));
G=simplify(CG-CG_0(0));

sCqp=size(Cqp);

for k=1:sCqp(1)
    [N,D]=numden(Cqp(k));
    N=simplifyFraction(N,'Expand',true)+0.7;
    terms=children(N);
    p=size(terms);
    for n=1:p(2)-1
        Ct(k,n)=terms(n)/D;
    end
end

sCt=size(Ct); %[xp,yp,thp,qp]
syms aux
C=zeros(4,4)*aux;
Qp=[x_op;theta_1p;theta_2p;theta_3p];

for n=1:sCt(1)
    kan=zeros(sCt(2),1);
    nn=1;
    for k=1:sCt(2)
        for m1=1:4
            [N0,D0]=numden(Ct(n,k));
            [N1,D1]=numden(Ct(n,k)/Qp(m1));
            eb=1;
            for kn=1:sCt(2)
                if kan(kn)==k
                    eb=0;
                end
            end
            if (D0==D1&&eb)
                kan(nn)=k;
                C(n,m1)=C(n,m1)+Ct(n,k)/Qp(m1);
                nn=nn+1;
            end
        end
    end
end
C;
G;
M;
Mp=M;
Cp=C;
Gp=G;
Ep=[1;0;0;0];
Fr=[-B_0,0,0,0;
    0,-B_1,0,0;
    0,0,-B_2,0;
    0,0,0,-B_3];
Cp=C-Fr;

Ax=simplify(-inv(Mp)*Cp);
Bx=simplify(inv(Mp)*Ep);
auxiliar=simplify(-inv(Mp)*Gp);
I=diag([1 1 1 1]);
Z=diag([0 0 0 0]);

A=[Ax,Z;I,Z];
B=[Bx;0;0;0;0];
Gfinal=[auxiliar;0;0;0;0];

SISTEMA=simplify(A*[x_op;theta_1p;theta_2p;theta_3p;x_o;theta1;theta2;theta3]+B*F+Gfinal);
S1=diff(SISTEMA(1),x_op)*x_op+diff(SISTEMA(1),theta_1p)*theta_1p+diff(SISTEMA(1),theta_2p)*theta_2p+diff(SISTEMA(1),theta_3p)*theta_3p+diff(SISTEMA(1),x_o)*x_o+diff(SISTEMA(1),theta1)*theta1+diff(SISTEMA(1),theta2)*theta2+diff(SISTEMA(1),theta3)*theta3;
S2=diff(SISTEMA(2),x_op)*x_op+diff(SISTEMA(2),theta_1p)*theta_1p+diff(SISTEMA(2),theta_2p)*theta_2p+diff(SISTEMA(2),theta_3p)*theta_3p+diff(SISTEMA(2),x_o)*x_o+diff(SISTEMA(2),theta1)*theta1+diff(SISTEMA(2),theta2)*theta2+diff(SISTEMA(2),theta3)*theta3;
S3=diff(SISTEMA(3),x_op)*x_op+diff(SISTEMA(3),theta_1p)*theta_1p+diff(SISTEMA(3),theta_2p)*theta_2p+diff(SISTEMA(3),theta_3p)*theta_3p+diff(SISTEMA(3),x_o)*x_o+diff(SISTEMA(3),theta1)*theta1+diff(SISTEMA(3),theta2)*theta2+diff(SISTEMA(3),theta3)*theta3;
S4=diff(SISTEMA(4),x_op)*x_op+diff(SISTEMA(4),theta_1p)*theta_1p+diff(SISTEMA(4),theta_2p)*theta_2p+diff(SISTEMA(4),theta_3p)*theta_3p+diff(SISTEMA(4),x_o)*x_o+diff(SISTEMA(4),theta1)*theta1+diff(SISTEMA(4),theta2)*theta2+diff(SISTEMA(4),theta3)*theta3;
S5=diff(SISTEMA(5),x_op)*x_op+diff(SISTEMA(5),theta_1p)*theta_1p+diff(SISTEMA(5),theta_2p)*theta_2p+diff(SISTEMA(5),theta_3p)*theta_3p+diff(SISTEMA(5),x_o)*x_o+diff(SISTEMA(5),theta1)*theta1+diff(SISTEMA(5),theta2)*theta2+diff(SISTEMA(5),theta3)*theta3;
S6=diff(SISTEMA(6),x_op)*x_op+diff(SISTEMA(6),theta_1p)*theta_1p+diff(SISTEMA(6),theta_2p)*theta_2p+diff(SISTEMA(6),theta_3p)*theta_3p+diff(SISTEMA(6),x_o)*x_o+diff(SISTEMA(6),theta1)*theta1+diff(SISTEMA(6),theta2)*theta2+diff(SISTEMA(6),theta3)*theta3;
S7=diff(SISTEMA(7),x_op)*x_op+diff(SISTEMA(7),theta_1p)*theta_1p+diff(SISTEMA(7),theta_2p)*theta_2p+diff(SISTEMA(7),theta_3p)*theta_3p+diff(SISTEMA(7),x_o)*x_o+diff(SISTEMA(7),theta1)*theta1+diff(SISTEMA(7),theta2)*theta2+diff(SISTEMA(7),theta3)*theta3;
S8=diff(SISTEMA(8),x_op)*x_op+diff(SISTEMA(8),theta_1p)*theta_1p+diff(SISTEMA(8),theta_2p)*theta_2p+diff(SISTEMA(8),theta_3p)*theta_3p+diff(SISTEMA(8),x_o)*x_o+diff(SISTEMA(8),theta1)*theta1+diff(SISTEMA(8),theta2)*theta2+diff(SISTEMA(8),theta3)*theta3;
S1_f=simplify(jacobian([S1,S2,S3,S4,S5,S6,S7,S8],[x_op, theta_1p,theta_2p,theta_3p,x_o,theta1,theta2,theta3]));
% 
U1=diff(SISTEMA(1),F)*F;
U2=diff(SISTEMA(2),F)*F;
U3=diff(SISTEMA(3),F)*F;
U4=diff(SISTEMA(4),F)*F;
U5=diff(SISTEMA(5),F)*F;
U6=diff(SISTEMA(6),F)*F;
U7=diff(SISTEMA(7),F)*F;
U8=diff(SISTEMA(8),F)*F;
U1_f=simplify(jacobian([U1,U2,U3,U4,U5,U6,U7,U8],[F]));

% 
x_o=0;
theta1=0;
theta2=0;
theta3=0;

theta_1p=0;
theta_2p=0;
theta_3p=0;
x_op=0;
F=0;

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

A_LINEAL=subs(S1_f)
B_LINEAL=subs(U1_f)

Co = ctrb(A_LINEAL,B_LINEAL);
unco = length(A_LINEAL) - rank(Co);
