function t = analyze(p_values, pearson_coefficients, alpha)
% Initialize counters
num_significant = 0;
num_non_significant = 0;
num_positive = 0;
num_negative = 0;

% Separate results for plotting
SignificantResults = [];
NonSignificantResults = [];

% Loop through data
for i = 1:length(p_values)
    if p_values(i) < alpha
        num_significant = num_significant + 1;
        if pearson_coefficients(i) > 0
            num_positive = num_positive + 1;
        elseif pearson_coefficients(i) < 0
            num_negative = num_negative + 1;
        end
        SignificantResults = [SignificantResults; pearson_coefficients(i), p_values(i)];
    else
        num_non_significant = num_non_significant + 1;
        NonSignificantResults = [NonSignificantResults; pearson_coefficients(i), p_values(i)];
    end
end

% Create the figure


% Plot significant results (if not empty)
if ~isempty(SignificantResults)
    semilogy(SignificantResults(:,1), SignificantResults(:,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    hold on;
end

% Plot non-significant results (if not empty)
if ~isempty(NonSignificantResults)
    semilogy(NonSignificantResults(:,1), NonSignificantResults(:,2), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

% Add a reference line for the significance threshold
yline(alpha, '--k', 'Significance Threshold (p=0.05)','FontWeight','bold','FontSize',10);

% Labels and title
xlabel('Pearson Correlation Coefficient','FontWeight','bold','FontSize',12);
ylabel('p-value (Log Scale)','FontWeight','bold','FontSize',12);
title('Pearson Correlation Coefficients and p-values','FontWeight','bold','FontSize',14);
ax = gca;
ax.FontSize = 10;           % bigger numbers
ax.FontWeight = 'bold';     % bold numbers

% Axis limits and grid
%     xlim([-1, 1]);  % Pearson correlation coefficient typically ranges from -1 to 1
%     ylim([1e-4, 1]);  % Adjust as needed for the range of p-values on the logarithmic scale
grid on;
box on;

% Add a legend
legend({
    sprintf('Significant (p < %.2f), N = %d (Positive = %d, Negative = %d)', alpha, num_significant, num_positive, num_negative), ...
    sprintf('Not Significant (p >= %.2f), N = %d', alpha, num_non_significant)
    }, 'Location', 'best','FontWeight','bold','FontSize',10);

hold off;
end
