clc; 
clear('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number_of_stations = 10;
lambda = 0.2;
a0 = 2*pi;
AR = 12;
alpha_geo = ((3.5)*pi)/180;
alpha_L0 = ((-1.5)*pi)/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all contributions? (Type (1) for yes or (2) for no %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotAll = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining Angle %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for station = 1:Number_of_stations
    theta(station, 1) = ((pi/2)/Number_of_stations)*station;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the trig functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odd = 1;
for column = 1:Number_of_stations
    for row = 1:Number_of_stations
    sin_functions(row, column) = sin(odd*theta(row, 1));
    end
    odd = odd+2;
end

for row = 1:Number_of_stations
    cos_functions(row, 1) = cos(theta(row, 1));
end
for row = 1:Number_of_stations
mu(row, 1) = ((a0)/(2*AR*(1+lambda))*(1+(lambda-1)*cos_functions(row, 1)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluating RHS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odd = 1;
for column = 1:Number_of_stations
    for row = 1:Number_of_stations
        A(row, column) = sin_functions(row, column)*(odd*mu(row, 1)+sin_functions(row, 1));
    end
    odd = odd+2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluating LHS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:Number_of_stations
    B(row, 1) = mu(row, 1)*(alpha_geo-alpha_L0)*sin_functions(row, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverting the matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = A\B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining CL % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cl = X(1,1)*pi*AR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wing Loading %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
column = 1;
odd = 1;
for row = 1:Number_of_stations
    for spanwise = -1:0.01:0
        Cl_normal(row,column) = X(row, 1)*sin(odd*acos(spanwise));
        Chord_percentage(1, column) = -(lambda-1)*(spanwise)+1;
        column = column + 1;
    end
    column = 1;
    odd = odd + 2;
end

for column = 1:101
    for row = 1:Number_of_stations
        if row == 1  
            Cl_normal_column_add(1, column) = Cl_normal(row, column);
        else
            Cl_normal_column_add(1, column) = Cl_normal_column_add(1, column) + Cl_normal(row, column);
        end
    end
end
    
for column = 1:101
        Cl_normal_plot(1,column) = ((2*(1+lambda)/(pi*X(1,1)))*(1/Chord_percentage(1,column)))*(Cl_normal_column_add(1,column));
end

for row = 1:Number_of_stations
    for column = 1:101
            Cl_normal_plot_comp(row,column) = ((2*(1+lambda)/(pi*X(1,1)))*(1/Chord_percentage(1,column)))*(Cl_normal(row,column));
    end
end 
spanwise = [-1:0.01:0];
plot(spanwise, Cl_normal_plot, 'k')
hold on
if plotAll == 1
    plot(spanwise, Cl_normal_plot_comp)
end
legend('cl/CL (y)')
title('Monoplane Wing Loading')
xlabel('y/(b/2)')
ylabel('cl/CL (y)')
