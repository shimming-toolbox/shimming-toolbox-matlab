function y_filt = filter_signal(probe_type, signal)
% Filter the input probe signal for drifts and high frequency glitches.
% 
% Input:
%   probe_type {'capacitive', 'pressure'}: Type of respiratory sensor
%   signal {str, vector}: Signal to filter. If signal is a string, it is 
%     considered as a file name. If it is a vector it is considered as 
%     the raw signal to filter.
% 

% TODO: implement feature based on probe_type
% signal = '/Users/julien/Desktop/capacitive_probe_baseline-capacitor';

%% Load signal
if isnumeric(signal)
    y = signal;

elseif ischar(signal)
    load(signal)
    y = data_cprobe;  % this is temporary: ultimately, "data_cprobe" should be called differently

end

%% Filter signal
x_first_sample_to_include_in_fitting = 300;
x_start_fitting = 600;

% timepoint = 1;
% count_zero_cross = 0;
% FS = stoploop({'Stop cprobe'}) ;
% while(~FS.Stop())
    % define signal from time 1 to time timepoint
%     signal_raw = data_cprobe(1:timepoint);

% x_full = 1:length(y);
x = 1:length(y);
x_sub = x_first_sample_to_include_in_fitting:length(y);

if length(x) > x_start_fitting
    start_fitting = true;
else
    start_fitting = false;
end

if start_fitting
    deg_poly = 3;
    [p, s, mu] = polyfit(x_sub, y(x_sub), deg_poly);
    % y_fit = (polyval(p, x_full) + mu(1)) * mu(2);
    y_fit = polyval(p, x_sub, [], mu);
    y_filt = y(x_sub) - y_fit;

else
    y_fit = y(x_sub);
    y_filt = y;
end

% <<< EXPERIMENTAL 
% Tried to fit with other functions:
% a = 1.8198153912055368E+04
% y_fit = a./sqrt(x_full)

% a = 1.0094966987333158E+00
% y_fit = power(a, x_full);
% >>>

do_plot = false;
if do_plot
    figure(10)
    % x = 1:timepoint;
    subplot(2,1,1)
    plot(x, y(x), 'k');
    hold on
    plot(x_sub, y_fit, 'r');
    hold off
    title('Raw signal.')  % Do we need that? --> Estimated frequency: ', num2str(meanfreq(y))])
    % legend('Raw signal', 'Fitted signal')
    ylabel('Frequency [KHz]')
    grid on
    if start_fitting
        subplot(2,1,2)
        plot(x_sub, y_filt, 'k');
        title(['Normalized signal (deg\_poly=' num2str(deg_poly) ')'])
        xlabel('Time [s]')  % TODO: fix label unit with sampling rate
        grid on
    end
end

% 
% function [zero_cross] = zero_cross(signal_in)
% product = signal_in(end) .* signal_in(end-1);
% if product > 0
%     zero_cross = 0;
% else
%     zero_cross = 1;
% end
% end
% 
% function [signal_hf, signal_lf] = separate_lf_hf(signal_in, x_start, x_stop, deg_poly)
% % Polynomial fitting for separating low and high frequency
% x = x_start:x_stop;
% p = polyfit(x, signal_in(x_start:x_stop), deg_poly);
% % define x variable based on signal_in
% x_full = 1:length(signal_in);
% signal_lf = polyval(p, x_full);
% signal_hf = signal_in - signal_lf;
% end
% 
% 
% function [x_zero, y_zero] = get_signal_passing_through_zero(signal_in)
% % this function will output a signal with half periods (assuming
% % quasi-sinusoinal shape). The half periods are detected by looking
% % at the passage through 0.
% for i = 2:length(signal_in)
%     % Search for values that pass thru 0 with ascending slope
%     if sign(signal_in(i) * signal_in(i - 1))
%         x_zero = i;
%         y_zero = signal_in(i);
%         break
%     end
% end
% end
% 
% 
% function deg_poly = update_deg_poly(deg_poly, timepoint)
% % Update polynomial degree based on number of timepoint. This is very empirical.
% if timepoint > 2000
%     deg_poly = 3;
% elseif timepoint > 1000
%     deg_poly = 2;
% elseif timepoint > 500
%     deg_poly = 1;
% elseif timepoint > 100
%     deg_poly = 0;
% end
% end
