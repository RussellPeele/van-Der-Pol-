%generate random fourier initial guess
%Newtons Method. How to decide initial guess? Runge kutta? 

x0 = .5;
y0 = 3;

% Period guess based on plot
T1 = 1000;
T2 = 2*pi
numPoints = 50;
mu = .1;
%w = T/(2*pi);

%Simulate transients:
[t,orbit] = ode45(@(t,y) odeRHS(t,y,mu),...
    linspace(0,T1,1000),[x0;y0]);

%Simulate attractor:
[t,orbit1] = ode45(@(t,y) odeRHS(t,y,mu),...
    linspace(0,T2,numPoints),orbit(end, :));
figure
hold on 
%plot(orbit(:,1), orbit(:,2))
plot(orbit1(:,1), orbit1(:,2))







%Compute Fourier coefficients from grid values:
B0x = ifft(orbit1(:,1)); 
B0y = ifft(orbit1(:,2));
K= numPoints;
%Coefficients are not in my prefered order - so fix it:
B00x = B0x(1:K/2+1);
B01x = B0x(K/2+2:end);
B00y = B0y(1:K/2+1);
B01y = B0y(K/2+2:end);
%C = [c_-n, c_-n+1,...c_-1, c_0,...c_n] as row vector 
Cx = [conj(B00x(end));B01x;B00x]';
Cy = [conj(B00y(end));B01y;B00y]';
goodguess=sum(Cx)

%order of series
NC = (size(Cx, 2)-1)/2;
azero = sum(Cx);

%van der pol Vector Field:
function dxdt = odeRHS(t,x,mu)
    dxdt = [0;0];
    dxdt(1) = x(2);
    dxdt(2) = mu*x(2)-mu*x(1)^2*x(2)-x(1);
    
end


