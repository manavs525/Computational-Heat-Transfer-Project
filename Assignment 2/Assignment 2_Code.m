clc
%Given values
A = 1; 
Ta = 1273;
sigma = 5.676e-8;
xspan = [0 0.2];
T0 = 290;
deltaT = 100;
errT = 0.0001;
k = 0;
%Initial guess for T2
T2 = 100;  
%convergence iteration 
while deltaT > errT
[x, T] = ode45(@heatT, xspan, T0, [], A,T2,Ta,sigma);
deltaT = abs(T2 - T(end));
T2 = T(end);
k = k+1;
end
%plot of temperature vs x
plot(x,T), xlabel('x(m)'), ylabel('T(K)');
fprintf("The value of T2 is : ");
disp(T2);
%function for heat transfer within a one-dimensional slab
function dT = heatT(x,T,A,T2,Ta,sigma)
dT = - sigma*(T2^4 - Ta^4)/(30*(1 + 0.002*T)*A);
end