% Given parameters
mu = 0.7;       % Mean
sigma = 0.1;    % Standard deviation

% Calculate Gamma parameters (k = shape, theta = scale)
k = (mu / sigma)^2;         % Shape parameter (α)
theta = (sigma^2) / mu;     % Scale parameter (β)

% Define x range (Gamma is defined for x > 0)
x = linspace(0, mu + 5*sigma, 1000);  % Cover up to mean + 5 std devs

% Compute PDF
pdf_values = gampdf(x, k, theta);

% Plot
figure;
plot(x, pdf_values, 'b-', 'LineWidth', 2);
title(sprintf('Gamma Distribution PDF (\\mu=%.2f, \\sigma=%.2f)', mu, sigma));
xlabel('x');
ylabel('Probability Density');
grid on;

% Annotate parameters
text(0.5*max(x), 0.9*max(pdf_values), ...
    sprintf('Shape (k) = %.2f\nScale (θ) = %.2f', k, theta), ...
    'FontSize', 10, 'VerticalAlignment', 'top');