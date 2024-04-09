function [positive_area, negative_area] = calculateAreaUnderComplexCurve(x, y)

% % Example vectors x and y
% x = linspace(0, 10, 1000); % Sample x values
% y = sin(x); % Sample y values (could be any curve)

% Find indices where y is positive and negative
positive_indices = y > 0;
negative_indices = y < 0;

% Split vectors x and y into positive and negative portions
x_positive = x(positive_indices);
y_positive = y(positive_indices);
x_negative = x(negative_indices);
y_negative = y(negative_indices);

% Calculate positive and negative areas using trapezoidal rule
positive_area = trapz(x_positive, y_positive);
negative_area = trapz(x_negative, y_negative);

% % Display results
% fprintf('Positive Area: %.2f\n', positive_area);
% fprintf('Negative Area: %.2f\n', negative_area);


end