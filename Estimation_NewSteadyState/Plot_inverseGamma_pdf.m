% Parameters
mu = 0.02;
sigma = 0.02;

% Calculate parameters
alpha = (2*sigma^2 + mu^2)/sigma^2;
beta = mu*(alpha - 1);

% Check if parameters are valid
if alpha <= 2
    error('For given mean and std dev, alpha must be > 2 for inverse gamma');
end

% Create x values (inverse gamma is defined for x > 0)
% Focus on small x values since mean is very small
x = linspace(0, 0.12, 1000);  % Adjust range as needed

% Calculate PDF
% Note: MATLAB doesn't have a built-in invgampdf, so we implement it
pdf_values = (beta^alpha)/gamma(alpha) * x.^(-alpha-1) .* exp(-beta./x);

% Plot
figure;
plot(x, pdf_values, 'LineWidth', 2);
title(sprintf('Inverse Gamma Distribution PDF (\\mu=%.2f, \\sigma=%.2f)', mu, sigma));
xlabel('x');
ylabel('Probability Density');
grid on;

% Display parameters
text(0.7*max(x), 0.9*max(pdf_values), ...
    sprintf('\\alpha = %.2f\n\\beta = %.4f', alpha, beta), ...
    'FontSize', 10, 'VerticalAlignment', 'top');

% Adjust x-axis if needed
xlim([0, 0.12]);  % Since mean is very small