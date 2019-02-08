%isolates and analyzes data from subset of points on a grid; takes output 
%from lymp_edges and analyzes one channel at a time

%New in v5: updated filtering algorithm 3 to identify outliers/ filter based on 95%
%CI and added visualization to check distribution of data vs. filtered data

%change vessel selection by boundary points or midpoints
%change outlier filtering by hard bound or number of sd's from mean or
%other metric

close all; clear all; clc;

%import matrix from csv - change filename/path as appropriate
coordinates = csvread('Vessel2-11.3-diameters.csv',1,0);



%%%%%%%%%%%%% VESSEL 1 %%%%%%%%%%%%%%


%restrict rectangle  of x and y coordinates to get vessel
%set boundaries
minx = 7;
maxx = 30;
miny = 20;
maxy = 100;


%filter the matrix around selected vessel
diameters = zeros(1,length(coordinates));
newcoordinates = zeros(length(coordinates),4);
for i=1:length(coordinates)
    if coordinates(i,4)>minx && coordinates(i,4)<maxx && coordinates(i,7)>miny && coordinates(i,7)<maxy
        %export vector of diameters
        diameters(i) = coordinates(i,1);
        %export vector of coordinates
        newcoordinates(i,1)=coordinates(i,2);
        newcoordinates(i,2)=coordinates(i,3);
        newcoordinates(i,3)=coordinates(i,5);
        newcoordinates(i,4)=coordinates(i,6);
    end
end

%save only the nonzero diameter values
diameters(diameters==0) = [];

%save only coordinate points corresponding to nonzero diameters
newcoordinates( ~any(newcoordinates,2), : ) = [];


%analyze distribution of diameters
%mean
diameters_mean = mean(diameters);
%standard deviation
diameters_stdev = std(diameters);
%standard error
diameters_sterr = diameters_stdev/sqrt(length(diameters)); 
tstat = tinv(0.975,length(diameters)-1); %95% CI with n-1 degrees of freedom for a sample of n diamter measurements
%threshold for refiltering
threshold_alg1 = 2*diameters_stdev; %2sd should capture 95% of the data, given large n to assume normality
threshold_alg2 = 0.01*diameters_mean;
threshold_alg3 = tstat * diameters_sterr;

%filtering algorithm 1
%define a new vector and filter outliers (+- 2sd from mean)
filtered_diameters = zeros(1,length(diameters));
filtered_coordinates = zeros(length(newcoordinates),4);
for i=1:length(diameters)
    if diameters(i) > diameters_mean-threshold_alg1 && diameters(i) < diameters_mean+threshold_alg1
        filtered_diameters(i)=diameters(i);
        filtered_coordinates(i,1)=newcoordinates(i,1);
        filtered_coordinates(i,2)=newcoordinates(i,2);
        filtered_coordinates(i,3)=newcoordinates(i,3);
        filtered_coordinates(i,4)=newcoordinates(i,4);
    end 
end

% %filtering algorithm 2
% %define a new vector and filter outliers based on the degree to which
% %including a measurement changes the mean of the sample (% difference)
% filtered_diameters = zeros(1,length(diameters));
% filtered_coordinates = zeros(length(newcoordinates),4);
% for i=1:length(diameters)
%     test_diameters = diameters([1:i-1,i+1:end]);
%     test_mean = mean(test_diameters);
%     test_std = std(test_diameters);
%     difference_std = test_std-diameters_stdev;
%     difference_mean = test_mean - diameters_mean;
%     
%     if abs(difference_mean)<threshold_alg2
%         filtered_diameters(i) = diameters(i);
%         filtered_coordinates(i,1)=newcoordinates(i,1);
%         filtered_coordinates(i,2)=newcoordinates(i,2);
%         filtered_coordinates(i,3)=newcoordinates(i,3);
%         filtered_coordinates(i,4)=newcoordinates(i,4);
%     end
%     
% 
% end


% %filtering algorithm 3
% %define a new vector and filter outliers (using a 95% confidence interval)
% filtered_diameters = zeros(1,length(diameters));
% filtered_coordinates = zeros(length(newcoordinates),4);
% for i=1:length(diameters)
%     if diameters(i) > diameters_mean-threshold_alg3 && diameters(i) < diameters_mean+threshold_alg3
%         filtered_diameters(i)=diameters(i);
%         filtered_coordinates(i,1)=newcoordinates(i,1);
%         filtered_coordinates(i,2)=newcoordinates(i,2);
%         filtered_coordinates(i,3)=newcoordinates(i,3);
%         filtered_coordinates(i,4)=newcoordinates(i,4);
%     end 
% end


%save only the new nonzero diameter values
filtered_diameters(filtered_diameters==0) = [];

%save only coordinate points corresponding to nonzero diameters
filtered_coordinates( ~any(filtered_coordinates,2), : ) = [];

%analyze new distribution of diameters
%mean
filtered_mean = mean(filtered_diameters);
%standard deviation
filtered_stdev = std(filtered_diameters);



%optional print statement (can comment this out if you don't need the mean
%and standard deviation displayed)
fprintf('Mean (filtered) =%s\n', filtered_mean);
fprintf('Standard deviation (filtered) =%s\n', filtered_stdev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plots filtered diameters

x=[filtered_coordinates(:,1) filtered_coordinates(:,2)];
y=[-filtered_coordinates(:,3) -filtered_coordinates(:,4)];
figure
plot(x',y', 'b')


%visualize distribution of diameters - another check to make sure filtering
%was reasonable

figure
histogram(diameters, 'BinWidth', 2) %original data distribution
xlim([0 max(diameters)+5])
title('Original Data Distribution')
figure
histogram(filtered_diameters, 'BinWidth', 2) %distribution of filtered data
xlim([0 max(diameters)+5])
title('Filterd Data Distribution')

%hold on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Uncomment below if you want to analyze or plot multiple vessels at once



% %%%%%%%%%%%%% VESSEL 2 %%%%%%%%%%%%%%
% 
% 
% %restrict rectangle  of x and y coordinates to get vessel
% %set boundaries
% minx2 = 125;
% maxx2 = 200;
% miny2 = 150;
% maxy2 = 200;
% 
% 
% %filter the matrix around selected vessel
% diameters2 = zeros(1,length(coordinates));
% newcoordinates2 = zeros(length(coordinates),4);
% for i=1:length(coordinates)
%     if coordinates(i,4)>minx2 && coordinates(i,4)<maxx2 && coordinates(i,7)>miny2 && coordinates(i,7)<maxy2
%         %export vector of diameters
%         diameters2(i) = coordinates(i,1);
%         %export vector of coordinates
%         newcoordinates2(i,1)=coordinates(i,2);
%         newcoordinates2(i,2)=coordinates(i,3);
%         newcoordinates2(i,3)=coordinates(i,5);
%         newcoordinates2(i,4)=coordinates(i,6);
%     end
% end
% 
% %save only the nonzero diameter values
% diameters2(diameters2==0) = [];
% 
% %save only coordinate points corresponding to nonzero diameters
% newcoordinates2( ~any(newcoordinates2,2), : ) = [];
% 
% 
% %analyze distribution of diameters
% %mean
% diameters_mean2 = mean(diameters2);
% %standard deviation
% diameters_stdev2 = std(diameters2);
% %threshold for refiltering
% threshold_alg1_2 = 2*diameters_stdev2;
% threshold_alg2_2 = 0.01*diameters_mean2;
% 
% %filtering algorithm 1
% % %define a new vector and filter outliers (+- 2sd from mean)
% % filtered_diameters2 = zeros(1,length(diameters2));
% % filtered_coordinates2 = zeros(length(newcoordinates2),4);
% % for i=1:length(diameters2)
% %     if diameters2(i) > diameters_mean2-threshold_alg1_2 && diameters2(i) <
% %     diameters_mean2+threshold_alg1_2
% %         filtered_diameters2(i)=diameters2(i);
% %         filtered_coordinates2(i,1)=newcoordinates2(i,1);
% %         filtered_coordinates2(i,2)=newcoordinates2(i,2);
% %         filtered_coordinates2(i,3)=newcoordinates2(i,3);
% %         filtered_coordinates2(i,4)=newcoordinates2(i,4);
% %     end 
% % end
% 
% %filtering algorithm 2
% %define a new vector and filter outliers based on the degree to which
% %including a measurement changes the mean of the sample (% difference)
% filtered_diameters2 = zeros(1,length(diameters2));
% filtered_coordinates2 = zeros(length(newcoordinates2),4);
% for i=1:length(diameters2)
%     test_diameters2 = diameters2([1:i-1,i+1:end]);
%     test_mean2 = mean(test_diameters2);
%     test_std2 = std(test_diameters2);
%     difference_std2 = test_std2-diameters_stdev2;
%     difference_mean2 = test_mean2 - diameters_mean2;
%     
%     if abs(difference_mean2)<threshold_alg2_2
%         filtered_diameters2(i) = diameters2(i);
%         filtered_coordinates2(i,1)=newcoordinates2(i,1);
%         filtered_coordinates2(i,2)=newcoordinates2(i,2);
%         filtered_coordinates2(i,3)=newcoordinates2(i,3);
%         filtered_coordinates2(i,4)=newcoordinates2(i,4);
%     end
%     
% 
% end
% 
% %save only the new nonzero diameter values
% filtered_diameters2(filtered_diameters2==0) = [];
% 
% %save only coordinate points corresponding to nonzero diameters
% filtered_coordinates2( ~any(filtered_coordinates2,2), : ) = [];
% 
% %analyze new distribution of diameters
% %mean
% filtered_mean2 = mean(filtered_diameters2);
% %standard deviation
% filtered_stdev2 = std(filtered_diameters2);
% 
% 
% 
% %optional print statement (can comment this out if you don't need the mean
% %and standard deviation displayed)
% fprintf('Vessel 2 Mean (filtered) =%s\n', filtered_mean2);
% fprintf('Vessel 2 Standard deviation (filtered) =%s\n', filtered_stdev2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %plots filtered diameters
% 
% x2=[filtered_coordinates2(:,1) filtered_coordinates2(:,2)];
% y2=[-filtered_coordinates2(:,3) -filtered_coordinates2(:,4)];
% plot(x2',y2', 'r')
% 
% hold on;
% 
% 
% 
% 
% 
% %%%%%%%%%%%%% VESSEL 3 %%%%%%%%%%%%%%
% 
% 
% %restrict rectangle  of x and y coordinates to get vessel
% %set boundaries
% minx3 = 126;
% maxx3 = 161;
% miny3 = 115;
% maxy3 = 150;
% 
% 
% %filter the matrix around selected vessel
% diameters3 = zeros(1,length(coordinates));
% newcoordinates3 = zeros(length(coordinates),4);
% for i=1:length(coordinates)
%     if coordinates(i,4)>minx3 && coordinates(i,4)<maxx3 && coordinates(i,7)>miny3 && coordinates(i,7)<maxy3
%         %export vector of diameters
%         diameters3(i) = coordinates(i,1);
%         %export vector of coordinates
%         newcoordinates3(i,1)=coordinates(i,2);
%         newcoordinates3(i,2)=coordinates(i,3);
%         newcoordinates3(i,3)=coordinates(i,5);
%         newcoordinates3(i,4)=coordinates(i,6);
%     end
% end
% 
% %save only the nonzero diameter values
% diameters3(diameters3==0) = [];
% 
% %save only coordinate points corresponding to nonzero diameters
% newcoordinates3( ~any(newcoordinates3,2), : ) = [];
% 
% 
% %analyze distribution of diameters
% %mean
% diameters_mean3 = mean(diameters3);
% %standard deviation
% diameters_stdev3 = std(diameters3);
% %threshold for refiltering
% threshold_alg1_3 = 2*diameters_stdev3;
% threshold_alg2_3 = 0.05*diameters_mean3;
% 
% % %filtering algorithm 1
% % %define a new vector and filter outliers (+- 2sd from mean)
% % filtered_diameters3 = zeros(1,length(diameters3));
% % filtered_coordinates3 = zeros(length(newcoordinates3),4);
% % for i=1:length(diameters3)
% %     if diameters3(i) > diameters_mean3-threshold_alg1_3 && diameters(i) < diameters_mean3+threshold_alg1_3
% %         filtered_diameters3(i)=diameters3(i);
% %         filtered_coordinates3(i,1)=newcoordinates3(i,1);
% %         filtered_coordinates3(i,2)=newcoordinates3(i,2);
% %         filtered_coordinates3(i,3)=newcoordinates3(i,3);
% %         filtered_coordinates3(i,4)=newcoordinates3(i,4);
% %     end 
% % end
% 
% %filtering algorithm 2
% %define a new vector and filter outliers based on the degree to which
% %including a measurement changes the mean of the sample (% difference)
% filtered_diameters3 = zeros(1,length(diameters3));
% filtered_coordinates3 = zeros(length(newcoordinates3),4);
% for i=1:length(diameters3)
%     test_diameters3 = diameters3([1:i-1,i+1:end]);
%     test_mean3 = mean(test_diameters3);
%     test_std3 = std(test_diameters3);
%     difference_std3 = test_std3-diameters_stdev3;
%     difference_mean3 = test_mean3 - diameters_mean3;
%     
%     if abs(difference_mean3)<threshold_alg2_3
%         filtered_diameters3(i) = diameters3(i);
%         filtered_coordinates3(i,1)=newcoordinates3(i,1);
%         filtered_coordinates3(i,2)=newcoordinates3(i,2);
%         filtered_coordinates3(i,3)=newcoordinates3(i,3);
%         filtered_coordinates3(i,4)=newcoordinates3(i,4);
%     end
%     
% 
% end
% 
% %save only the new nonzero diameter values
% filtered_diameters3(filtered_diameters3==0) = [];
% 
% %save only coordinate points corresponding to nonzero diameters
% filtered_coordinates3( ~any(filtered_coordinates3,2), : ) = [];
% 
% %analyze new distribution of diameters
% %mean
% filtered_mean3 = mean(filtered_diameters3);
% %standard deviation
% filtered_stdev3 = std(filtered_diameters3);
% 
% 
% 
% %optional print statement (can comment this out if you don't need the mean
% %and standard deviation displayed)
% fprintf('Vessel 3 Mean (filtered) =%s\n', filtered_mean3);
% fprintf('Vessel 3 Standard deviation (filtered) =%s\n', filtered_stdev3);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %plots filtered diameters
% 
% x3=[filtered_coordinates3(:,1) filtered_coordinates3(:,2)];
% y3=[-filtered_coordinates3(:,3) -filtered_coordinates3(:,4)];
% plot(x3',y3', 'g')
% 
% hold off
