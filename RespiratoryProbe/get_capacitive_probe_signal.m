function get_capacitive_probe_signal()
% Main function that retrieves capacitive probe signal in real time and
% filter it appropriately.
% 
%   FS: Structure that gathers info from probe
%   data_cprobe: signal from probe

load('../test/capacitive_probe_phantom.mat');  % for debugging
n_samples_before_start_fitting = 50;
deg_poly = 0;

timepoint = 1;
FS = stoploop({'Stop cprobe'}) ;
while(~FS.Stop())
    % define signal from time 1 to time timepoint
    signal_norm = data_cprobe(1:timepoint);

    % Update polynomial degree for normalization
    deg_poly = update_deg_poly(deg_poly, timepoint);
    
    % Display stuff
    if timepoint > n_samples_before_start_fitting
        [signal_norm, signal_fit] = update_fit(signal_norm, 1, deg_poly);
        % Display
        figure(1)
        x = 1:timepoint;
        subplot(2,1,1)
        plot(x, data_cprobe(x), 'k');
        hold on
        plot(x, signal_fit, 'r');
        hold off
        title('Raw signal')
        ylabel('Frequency [KHz]')
        grid on
        subplot(2,1,2)
        plot(x, signal_norm, 'k');
        title(['Normalized signal (deg\_poly=' num2str(deg_poly) ')'])
        xlabel('Time [s]')
        grid on
    else
        % Display
        figure(1)
        x = 1:timepoint;
        subplot(2,1,1)
        plot(x, data_cprobe(x), 'k');
        grid on
    end
    timepoint = timepoint + 1;
 
end
end


function [signal_norm, signal_fit] = update_fit(signal_in, x_start, deg_poly)
% Update the fit and display fig
x = x_start:length(signal_in);
p = polyfit(x, signal_in(x_start:end), deg_poly);
% define x variable based on signal_in
x_full = 1:length(signal_in);
signal_fit = polyval(p, x_full);
signal_norm = signal_in - signal_fit;
end


function [x_zero, y_zero] = get_signal_passing_through_zero(signal_in)
% this function will output a signal with half periods (assuming
% quasi-sinusoinal shape). The half periods are detected by looking
% at the passage through 0.
for i = 2:length(signal_in)
    % Search for values that pass thru 0 with ascending slope
    if sign(signal_in(i) * signal_in(i - 1))
        x_zero = i;
        y_zero = signal_in(i);
        break
    end
end
end


function deg_poly = update_deg_poly(deg_poly, timepoint)
% Update polynomial degree based on number of timepoint. This is very empirical.
if timepoint > 2000
    deg_poly = 3;
elseif timepoint > 1000
    deg_poly = 2;
elseif timepoint > 500
    deg_poly = 1;
elseif timepoint > 100
    deg_poly = 0;
end
end