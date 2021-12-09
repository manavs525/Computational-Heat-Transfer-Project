clc;
t = [0:5.5:121];                %range of t(121 gives an integer value and is close to 120)
x = [0:0.05:0.85];              %range of distance in x

%Explicit Method
m = length(t);
n = length(x);
T = zeros(m,n);
%Initial Temperature
for i=1:m
    for j=1:n
        T(i,j)=20;
    end
end
%Calculating Temperature at surface and interior nodes
for i=2:m
    T(i,1)=(24.94+T(i-1,2))/2+0.5*T(i-1,1);                 %Equation for surface node(Calculations and derivation are shown in report)
    for j=2:n-1
       T(i,j) = (T(i-1,j+1) + T(i-1,j-1))/4+0.5*T(i-1,j);   %Equation for interior nodes
    end
end

disp("THE EXPLICIT METHOD");
disp("The surface temperature after 2 min is :");
disp(T(i,1));
disp("The temperature at x=150mm after 2 min is :");
disp(T(i,4));
disp(T);                            %Matrix output to get idea for implicit method

figure(1)
surf(x,t,T);
colorbar

%IMPLICIT Method
A = zeros(16,16);                   %Coefficient matrices
B = zeros(16,1);
for i = 1:16
    for j = 1:16
        if(i == j)
            A(i,j) = 6;
        elseif(j == i-1 || j == i+1)
                A(i,j) = -1;
        end
    end
end
A(1,1)=3;
T1 = zeros(23,16);
for i=1:23                              %Initial Temperature
    for j = 1:16
    T1(i,j)=20;
    end
end

for i = 2:23                            %Calculation of B Matrix and Temperature matrix
    for j = 1:15
        B(j,1)=4*T1(i-1,j);
    end
    B(1,1)= 25.68 + 2*T1(i-1,1);
    B(16,1) = 4*T1(i-1,16) + 20;
       
    T1(i,:) = linsolve(A,B);
end
disp("THE IMPLICIT METHOD");
disp("The surface temperature after 2 min is :");
disp(T1(i,1));
disp("The temperature at x=150mm after 2 min is :");
disp(T1(i,4));

figure(2)
t = [0:5.5:121];
x = [0:0.05:0.75];
surf(x,t,T1);
colorbar
