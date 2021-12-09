clc
t = [0:1:15];
T = zeros(45,1);
T1 = 293*ones(45,1);
%Unsteady Heat transfer
[t, T] = ode45(@(t,T)myode(t,T), t, T1);
tempgrid = steadytemperaturedistribution(T(16,:));

%Plot
X = [0:1/9:1];
Y = [0:0.05:0.2];
s1 = surf(X,Y,tempgrid);
s1.FaceColor = 'texturemap';
colorbar;
xlabel("X");
ylabel("Y");
zlabel("Temperature");
hold on;

%Steady state heat transfer
[A1,B1] = Calc();
Tempvec = linsolve(A1,B1);
TEMPGRID = steadytemperaturedistribution(Tempvec);
s2 = surf(X,Y,TEMPGRID);
s2.FaceColor = 'interp';
colorbar;
hold off;

%function for unsteady heat conduction
function F = myode(t,T)
    [A,B] = Calc();
    F = A*T-B;
end

%function to calculate A and B coeffecient matrix
function [A,B] = Calc()
    dx = 1/9;
    dy = 0.05;
    k = 401;
    h = 10;
    A = zeros(45,45);
    B = zeros(45,1);

    %coeffecients for top node(node 1)
    A(1,1)= -h*dx - k*dy/dx - k*dx/dy;
    A(1,2)= k*dy/(2*dx);
    A(1,10)= k*dx/dy;
    B(1)= -303*k*dy/(2*dx) - h*dx*293;

    %coefficients for top nodes(nodes 2-8)
    for i=2:8
        A(i,i)= -h*dx - k*dy/dx - k*dx/dy;
        A(i,i+1)= k*dy/(2*dx);
        A(i,i-1)= k*dy/(2*dx);
        A(i,i+9)= k*dx/dy;
        B(i)= -h*dx*293;
    end

    %coefficients for top right corner node(nodes 9)
    A(9,9)= (-h*dx - k*dy/dx - k*dx/dy)/2;
    A(9,8)= k*dy/(2*dx);
    A(9,18)= k*dx/(2*dy);
    B(9)= -h*dx*293/2;

    %coefficients for interior nodes and right nodes(nodes 10-36)
    for i=10:36
        if i==10||i==19||i==28
            A(i,i)= -2/dx^2 - 2/dy^2;
            A(i,i+1)= 1/dx^2;
            A(i,i+9)= 1/dy^2;
            A(i,i-9)= 1/dy^2;
            B(i)= -1/dx^2*303;
        elseif i==18||i==28||i==36
            A(i,i)= -2/dx^2 - 2/dy^2;
            A(i,i-1)= 2*1/dx^2;
            A(i,i+9)= 1/dy^2;
            A(i,i-9)= 1/dy^2;
            B(i)=0;
        else
            A(i,i)= -2/dx^2 - 2/dy^2;
            A(i,i+1)= 1/dx^2;
            A(i,i-1)= 1/dx^2;
            A(i,i+9)= 1/dy^2;
            A(i,i-9)= 1/dy^2;
            B(i)=0;
        end
    end

    %coefficients for bottom node(node 37)
    A(37,37)= -h*dx - k*dy/dx - k*dx/dy;
    A(37,38)= k*dy/(2*dx);
    A(37,28)= k*dx/dy;
    B(37)= -303*k*dy/(2*dx) - h*dx*293;

    %coefficients for bottom nodes(node 38-44)
    for i=38:44
        A(i,i) = -h*dx - k*dy/dx - k*dx/dy;
        A(i,i+1) = k*dy/(2*dx);
        A(i,i-1) = k*dy/(2*dx);
        A(i,i-9) = k*dx/dy;
        B(i) = -h*dx*293;
    end

    %coefficients for bottom right corner node(node 45)
    A(45,45) = (-h*dx - k*dy/dx - k*dx/dy)/2;
    A(45,44) = k*dy/(2*dx);
    A(45,36) = k*dx/(2*dy);
    B(45)= -h*dx*293/2;
end

%function to map temperature vector into a 2D grid
function Tempgrid = steadytemperaturedistribution(Tempvec)
    Tempgrid = zeros(5,10);
    k1 = 1;
    for i1 = 1:5
        for j = 1:10
            if j==1
                Tempgrid(i1,j) = 303;
            else
                Tempgrid(i1,j) = Tempvec(k1);
                k1 = k1+1;
            end
        end
    end
end
