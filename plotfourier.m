%plot fourier
theta = linspace(-pi, pi, numPoints);
xarray = zeros(1,numPoints);
yarray = zeros(1,numPoints);
for k=1:numPoints
    xarray(k)= evaluateCFourier(theta(k),order);
    yarray(k)= evaluateCFourier(theta(k),order);
end
figure 
hold on 
plot(xarray, yarray)