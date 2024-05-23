clc; clear all; close all;

syms m M L g gam real
syms q1 q2 real
syms dq1 dq2 dq12 dq22 real
 
q=[q1;q2];
dq=[dq1;dq2];

% Kinetic Energy
rM=L*[-sin(q1);cos(q1)];
rm=L*[-sin(q1);cos(q1)]+L*[-sin(q1-q2);-cos(q1-q2)];
VM=L*dq1*[cos(q1);sin(q1)];
Vm=L*dq1*[cos(q1);sin(q1)]+L*(dq1-dq2)*[cos(q1-q2);-sin(q1-q2)];
T= 1/2*M*VM.'*VM+1/2*m*Vm.'*Vm;

% Potential Energy
V=M*g*L*cos(q1-gam)+m*g*(L*cos(q1-gam)-L*cos(q1-gam-q2));

L=T-V;

dL_ddq=jacobian(L,dq).';

Mass_Matrix= jacobian(dL_ddq,dq);
N = jacobian(dL_ddq , q);
B = N*dq - jacobian(L,q).';

Mass_Matrix = simplify(Mass_Matrix);
B = simplify(B);

%angular momentum
MP2_1=M*cross(L*[-sin(-q1);cos(-q1);0],[VM;0])+m*cross(L*[-sin(-q1);cos(-q1);0]+L*[-sin(-q1+q2);-cos(-q1+q2);0],[VM;0]);
MH_1=m*cross(L*[-sin(-q1+q2);-cos(-q1+q2);0],[VM;0]);
MP2_2=M*cross(L*[-sin(-q1);cos(-q1);0],L*dq12*[cos(-q1);sin(-q1);0])+m*cross(L*[-sin(-q1);cos(-q1);0]+L*[-sin(-q1+q2);-cos(-q1+q2);0],L*dq12*[cos(-q1);sin(-q1);0]+L*(dq12-dq22)*[cos(-q1+q2);-sin(-q1+q2);0]);
MH_2=m*cross(L*[-sin(-q1+q2);-cos(-q1+q2);0],L*dq12*[cos(-q1);sin(-q1);0]+L*(dq12-dq22)*[cos(-q1+q2);-sin(-q1+q2);0]);

F=[MP2_2(3)-MP2_1(3);MH_2(3)-MH_1(3)];
AA=jacobian(F,[dq12 dq22]);
BB=AA*[dq12;dq22]-F;
AA=simplify(AA)
BB=simplify(BB)