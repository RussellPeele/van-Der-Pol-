%testing differential suboperators for f(x,y) = (y, mu*y -mu*x^2y-x)

%phase condition x(0)=0 




N=length(Cx);


%Newtons Method. How to decide initial guess? Runge kutta? 
azero = 1.260122998207854;
bzero = 1.632578584872354;
a = Cx;
b = Cy;

order = (length(a)-1)/2;

w = T/(2*pi);
DF=buildDF(w,a,b,mu);
L=w*makeML(N);
epsilon = 10^-6;
defect = norm(evalF(w,a, b,azero,mu,N), "inf");





 xguess=[w,a,b];

 %test DF 

 

for j=1:r
    deltabar(j)=delta;
end

while norm(defect, "inf") > epsilon
    h=xguess(1);
    u=xguess(2:(r-1)/2+1);
    v=xguess((r-1)/2+2:end);
    Fwab = evalF(h,u,v,azero, mu,N)';
    defect = norm(Fwab, "inf");
   
    correction = transpose(-DF\Fwab); 
    xguess = xguess + correction; 
end





