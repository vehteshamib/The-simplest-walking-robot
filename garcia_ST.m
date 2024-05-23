function dz = garcia_ST(t,z,t_init)

q1=z(1);  dq1=z(2);  
q2=z(3);  dq2=z(4);  
global m M L g gam

AA =[ L^2*(M + 2*m + 2*m*cos(2*q1 - q2)), -L^2*m*(cos(2*q1 - q2) + 1)
	                  -L^2*m*(cos(2*q1 - q2) + 1),                       L^2*m];
B =[L*M*g*sin(gam - q1) + L*g*m*sin(gam - q1) - 2*L^2*dq1^2*m*sin(2*q1 - q2) - L^2*dq2^2*m*sin(2*q1 - q2) - L*g*m*sin(gam - q1 + q2) + 2*L^2*dq1*dq2*m*sin(2*q1 - q2)
                                                                                                                m*sin(2*q1 - q2)*L^2*dq1^2 + g*m*sin(gam - q1 + q2)*L];

ddq=-inv(AA)*B;
dz=[dq1 ddq(1) dq2 ddq(2)].';
