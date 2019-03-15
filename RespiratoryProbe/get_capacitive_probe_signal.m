function get_capacitive_probe_signal()
% Main function that retrieves capacitive probe signal in real time and
% filter it appropriately.
% List of functions dealing with the capacitive probe.

% uigetfile;

clc; clear all; close all;
load('../test/capacitive_probe_phantom.mat');
n_samples_before_start_fitting = 300;

deg_poly = 3;
preliminary_fit = 1;

timepoint = 1;
FS = stoploop({'Stop cprobe'}) ;
while(~FS.Stop())

    % define signal from time 1 to time timepoint
    filtered_signal = data_cprobe(1:timepoint);

    % Wait a number of samples before fitting
    if timepoint > n_samples_before_start_fitting

        % preliminary fit before we know the number of periods
        if preliminary_fit
            filtered_signal = update_fit(filtered_signal, 1, 1)
        end

        % Search for values that pass thru 0 with ascending slope
        if sign(filtered_signal(timepoint) * filtered_signal(timepoint - 1)) == -1 && filtered_signal(timepoint) > filtered_signal(timepoint - 1)

            % Get signal with complete periods
            x_start = get_complete_periods(filtered_signal);

            % Update fit
            filtered_signal = update_fit(filtered_signal, x_start, deg_poly);

            % stop preliminary fit, because now we will only fit at
            % complete periods
            preliminary_fit = 0
        end
    end
    timepoint = timepoint + 1;
end
end

function signal_out = update_fit(signal_in, x_start, deg_poly)
    x = x_start:length(signal_in);
    p = polyfit(x, signal_in(x_start:end), deg_poly);
    % define x variable based on signal_in
    x_full = 1:length(signal_in);
    y = polyval(p, x_full);
    signal_out = signal_in - y;

    % Display stuff
    figure(1)
    subplot(2,1,1)
    plot(x_full, y, 'r');
    subplot(2,1,2)
    plot(x_full, signal_out, 'k');
    title('signal - polyfit 4')
    xlabel('Time [s]')
    ylabel('Frequency [KHz]')
    grid on    
end


function x_start = get_complete_periods(signal_in)
    % this function will output a signal with complete periods (assuming
    % quasi-sinusoinal shape). The number of periods is detected by looking
    % at the passage through 0.
    for i = 2:length(signal_in)
        % Search for values that pass thru 0 with ascending slope
        if sign(signal_in(i) * signal_in(i - 1)) == -1 && signal_in(i) > signal_in(i - 1)
            x_start = i;
            break
        end
    end
end

