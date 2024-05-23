function [gstop, isterminal,direction]=collision(t,z)
q1 = z(1); q2 = z(3); 
gstop =(2*q1-q2);

if (q2>-0.05) %no collision detection for foot scuffing
    isterminal = 0;
else
    isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
end
direction=-1;