% Inputs to play with
% Initialize cursors
% Note behaviour when sum(abs(post_cursors)) > abs(main_cursor) - sum(abs(pre_cursors))
post_cursors = [-1.1, 2.8, -4.1, 3.3, 2.1, -13.6, 13.3, -20.5];
main_cursor = 100;
pre_cursors = [30, -10, 10];
% Increase decimation_factor if you run out of memory
decimation_factor = 1e3;
% Increase sequence_length if equalizer constants need more transmitions to converge
sequence_length = 1e8;
% Increase step_size to make equalizer converge faster, too much will ruin convergence
step_size = 1e-6;



decimated_idx = 1;
% Generate a random sequence of 1's and -1's
bits = 2 * randi([0, 1], 1, sequence_length) - 1;

% Pad with zeros
padding_length = max(length(pre_cursors), length(post_cursors));
bits = [zeros(1, padding_length), bits, zeros(1, padding_length)];

% Generate receiver noisy voltages using the cursors
voltages = generate_voltages(bits, pre_cursors, post_cursors, main_cursor);

% Initialize variables for DFE
estimated_voltages = voltages;
estimated_bits = zeros(1, length(bits));
coeffs = zeros(1, length(post_cursors));

% Initialize storage for coefficients history
coeffs_history = zeros(length(post_cursors), floor(sequence_length / decimation_factor));

figure;
xlabel('time [bit sample interval]');
ylabel('Coefficient Value');
title('DFE Coefficients Over Time');
hold on;

% Analog-Sign DFE processing
start_idx = padding_length + 1;
end_idx = length(estimated_voltages) - padding_length;
for n = start_idx:end_idx
    for k = 1:length(coeffs)
        estimated_voltages(n) = estimated_voltages(n) - coeffs(k) * estimated_bits(n-k);
    end
    estimated_bits(n) = sign(estimated_voltages(n));
    % Update coefficients
    if estimated_bits(n) == 1
        estimated_error_signal = (estimated_voltages(n) - main_cursor);
    else
        estimated_error_signal = (estimated_voltages(n) + main_cursor);
    end
    for k = 1:length(coeffs)
        coeffs(k) = coeffs(k) + step_size * estimated_error_signal * estimated_bits(n - k);
    end
    
    % Store the updated coefficients for plotting
    if mod(n, decimation_factor) == 0
        coeffs_history(:, decimated_idx) = coeffs';
        decimated_idx = decimated_idx + 1;
    end
    
    % Print progress every million estimations
    if mod(n, 1e6) == 0
        fprintf('Processed %d million bits.\n', n / 1e6);
    end
    
    % Update plot every 10,000 estimations (or choose another interval)
    if mod(n, 1e6) == 0  % Reduced for demonstration
    cla; % Clear previous plots
    for k = 1:length(coeffs)
        plot((1:decimated_idx-1) * decimation_factor, coeffs_history(k, 1:decimated_idx - 1), 'DisplayName', ['Coeff ' num2str(k)]);
    end
    legend; % Show legend

    % Determine plotting range to fit decimated data
    min_x = 0;
    max_x = (decimated_idx - 1) * decimation_factor;

    % Plot horizontal lines for post-cursors, but restrict their range
    for k = 1:length(post_cursors)
        line([min_x, max_x], [post_cursors(k), post_cursors(k)], 'Color', 'k', 'LineStyle', '--', 'DisplayName', ['Post-cursor ' num2str(k)]);
    end

    drawnow;  % Redraw the plot with the updated coefficients
end
end
%drawnow;
% Calculate Errors Per Million (EPM)
errors = sum(bits(start_idx:end_idx) ~= estimated_bits(start_idx:end_idx));
EPM = (errors / sequence_length) * 1e6;
disp(['Errors: ', num2str(errors)]);
disp(['Errors Per Million (EPM): ', num2str(EPM)]);


% Generates noisy voltages that the receiver would measure when a transmitter tries to send the bits from function input
function voltages = generate_voltages(bits, pre_cursors, post_cursors, main_cursor)
    % Initialize voltages to zero
    voltages = zeros(1, length(bits));

    % Compute voltages using pre-cursors, main cursor, and post-cursors
    padding_length = max(length(pre_cursors), length(post_cursors));
    start_idx = padding_length + 1;
    end_idx = length(voltages) - padding_length;
    
    % Loop over the none pad part
    for n = start_idx:end_idx
        % Add the effect of the main cursor
        voltages(n) = voltages(n) + bits(n) * main_cursor;

        % Add the effect of pre-cursors
        for k = 1:length(pre_cursors)
            voltages(n) = voltages(n) + bits(n + k) * pre_cursors(k);
        end

        % Add the effect of post-cursors
        for k = 1:length(post_cursors)
            voltages(n) = voltages(n) + bits(n - k) * post_cursors(k);
        end
    end
end
