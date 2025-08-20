clear;
close all;

% Parameters
mu = 0.5;
sigma = 0.1;

% Calculate alpha and beta
v = sigma^2;
alpha = ((1 - mu) / v - 1 / mu) * mu^2;
beta = alpha * (1 / mu - 1);

% Create x values
x = linspace(0, 1, 1000);

% Calculate PDF
pdf_values = betapdf(x, alpha, beta);

% Plot
figure;
plot(x, pdf_values, 'LineWidth', 2);
title(sprintf('Beta Distribution PDF (\\mu=%.2f, \\sigma=%.2f)', mu, sigma));
xlabel('x');
ylabel('Probability Density');
grid on;

% Display parameters
text(0.1, max(pdf_values)*0.9, sprintf('\\alpha = %.2f\n\\beta = %.2f', alpha, beta));