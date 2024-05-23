clear all
close all
clc
format long
global m M L g gam
m=0.0025; M=1; L=1; g=10; gam=0.;
z_init=[0.2 -0.2 0.4 0].';
t_init=0;
 
for step = 1:1
   
step
    options = odeset('abstol',1e-9,'reltol',1e-9,'events',@collision);
    dif=@(t,z) garcia_ST(t,z,t_init);
    [t,z] = ode45(dif,[t_init 10],z_init,options);
        
    z_init=garcia_HS(z_init);
    
    plot(t,z(:,1),'b')
    hold on
    plot(t,z(:,3),'r')      
        
end
