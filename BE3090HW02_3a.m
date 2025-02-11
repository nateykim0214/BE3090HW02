% Parameters for the first normal distribution
mu1 = 5; % mean
sigma1 = 1; % standard deviation

% Parameters for the second normal distribution
mu2 = 7; % mean
sigma2 = 1; % standard deviation

% Define a range of x values for the distributions
x = 0:.01:10; % From 1 to 11 (to cover both distributions)

% Calculate the probability density function (PDF) for both distributions
pdf1 = normpdf(x, mu1, sigma1);
pdf2 = normpdf(x, mu2, sigma2);

% Combine the two PDFs
combined_pdf = pdf1 + pdf2;

% Normalize the combined PDF
normalized_pdf = combined_pdf / trapz(x, combined_pdf); % trapz for numerical integration

% Plot the individual distributions and the normalized combined PDF
figure;
plot(x, pdf1, 'b-', 'LineWidth', 2); % Distribution 1
hold on;
plot(x, pdf2, 'r-', 'LineWidth', 2); % Distribution 2
hold on;
plot(x, normalized_pdf, 'g-', 'LineWidth', 2); % Normalized combined PDF

% Labels and legend
xlabel('x');
ylabel('Probability Density');
title('Normal Distributions and Normalized Combined PDF');
legend('Distribution 1 (mu = 5)', 'Distribution 2 (mu = 7)', 'Normalized Combined PDF');
grid on;
