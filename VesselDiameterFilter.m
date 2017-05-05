%isolates and analyzes data from subset of points on a grid; takes output 
%from lymp_edges and analyzes one channel at a time


close all; clear all; clc;

%import matrix from csv - change filename/path as appropriate
coordinates = csvread('test_lymph2_AMTtrial-diameters.csv',1,0);


%restrict rectangle  of x and y coordinates to get vessel
%set boundaries
minx = 25;
maxx = 75;
miny = 75;
maxy = 125;


%filter the matrix
diameters = zeros(1,length(coordinates));
for i=1:length(coordinates)
    if coordinates(i,2)>minx && coordinates(i,2)<maxx && coordinates(i,3)>miny && coordinates(i,3)<maxy
        %export vector of diameters
        diameters(i) = coordinates(i,1);
    end
end

%save only the nonzero diameter values
diameters(diameters==0) = [];

diameters

%analyze distribution of diameters
%mean
mean = mean(diameters)
%standard deviation
stdev = std(diameters)

%optional print statement (can comment this out if you don't need the mean
%and standard deviation displayed)
fprintf('Mean =%s\n', mean);
fprintf('Standard deviation =%s\n', stdev);

