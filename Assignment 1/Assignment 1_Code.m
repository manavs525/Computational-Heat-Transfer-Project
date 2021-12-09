% L = 0.2, b = 0.04, dx = 0.0105, T(0) = 100, ;
% Final Equation used is
%[1-(m-0.5)*0.0525]T(m-1)+[-2+0.105*m-0.00046]T(m)+[1-(m+0.5)*0.0525]T(m+1)= -0.0115
clc
temp = zeros(19,1);
temp(1,1) = 100;
A = zeros(20);
B = zeros(20,1);
A(1,1) = -2 + 0.105 - 0.00046;
A(1,2) = 1 - 1.5*0.0525;
B(1,1) = -97.375;
for m = 2:19
    A(m,m-1)= 1 - (m - 0.5)*0.0525;
    A(m,m) = -2 + 0.105*m - 0.00046;
    A(m,m+1) = 1 - (m + 0.5)*0.0525;
    B(m,1) = -0.0115;
end
%The final coefficients are obtained by writing an energy balance on the volume element of length dx/2
A(20,19) = 1;
A(20,20) = -1.009;
B(20,1) = -0.189;
tnodes = linsolve(A,B);
for i=2:20
    temp(i) = tnodes(i-1);
end
disp("Temperature Values at the Nodes: ");
disp(temp);

x = (0:0.0105:0.2);
plot(x,temp)
title('Plot of Temperature v/s Distance');
xlabel('Distance');
ylabel('Temperature');