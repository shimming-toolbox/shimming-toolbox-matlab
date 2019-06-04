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
count_zero_cross = 0;
FS = stoploop({'Stop cprobe'}) ;
while(~FS.Stop())
    % define signal from time 1 to time timepoint
    signal_raw = data_cprobe(1:timepoint);

    % Update polynomial degree for normalization
    deg_poly = update_deg_poly(deg_poly, timepoint);
    
    % Display stuff
    if timepoint > n_samples_before_start_fitting   
        [signal_hf, signal_lf] = separate_lf_hf(signal_raw, 1, stop_time, deg_poly);
        if zero_cross(signal_hf)
            count_zero_cross = count_zero_cross + 1;
            if rem(count_zero_cross, 2) == 0
                stop_time = timepoint;
            end
        end
        % Display
        figure(1)
        x = 1:timepoint;
        subplot(2,1,1)
        plot(x, data_cprobe(x), 'k');
        hold on
        plot(x, signal_lf, 'r');
        hold off
        title(['Raw signal estimated frequency is : ', num2str(meanfreq(signal_raw))])
        ylabel('Frequency [KHz]')
        grid on
        subplot(2,1,2)
        plot(x, signal_hf, 'k');
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
        stop_time = timepoint;
    end
    timepoint = timepoint + 1;
 
end
end

function [zero_cross] = zero_cross(signal_in)
product = signal_in(end) .* signal_in(end-1);
if product > 0
    zero_cross = 0;
else
    zero_cross = 1;
end
end

function [signal_hf, signal_lf] = separate_lf_hf(signal_in, x_start, x_stop, deg_poly)
% Polynomial fitting for separating low and high frequency
x = x_start:x_stop;
p = polyfit(x, signal_in(x_start:x_stop), deg_poly);
% define x variable based on signal_in
x_full = 1:length(signal_in);
signal_lf = polyval(p, x_full);
signal_hf = signal_in - signal_lf;
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