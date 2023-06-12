% Extract Data
data = Extract_Half_Car_Rig_Data();

Fs = 5000;  % sampling frequency

% Create tiled layout
tiledlayout(8,4)

% Iterate over all data sets
for i=1:size(data,2)
    % Compute FFT of the data
    fft_data = fft(data(i).rawdata(:,4));
    L = length(data(i).rawdata(:,4));  % length of the signal
    f = Fs*(0:(L/2))/L;  % frequency vector
    P2 = abs(fft_data/L);  % double-sided spectrum
    P1 = P2(1:L/2+1);  % single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);

    % Plot the single-sided amplitude spectrum in the next tile
    nexttile;
    semilogy(f,P1)
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
    xlabel('Frequency (f)')
    ylabel('|P1(f)| in dB')
    
    % Set the x-axis limit to 100
    xlim([0 20])
end




%% Extract Data
data = Extract_Half_Car_Rig_Data()

fft_data = fft(data(29).rawdata(:,4));
Fs = 5000;  % sampling frequency
L = 20000;  % length of the signal
f = Fs*(0:(L/2))/L;  % frequency vector
P2 = abs(fft_data/L);  % double-sided spectrum
P1 = P2(1:L/2+1);  % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

% Now, to plot the single-sided amplitude spectrum
figure;
semilogy(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Frequency (f)')
ylabel('|P1(f)| in dB')

%% Plot Data - fft, tiled layout 
tiledlayout(8,4)
for i=1:size(data,2);
    nexttile;
    plot(data(i).catagoriseddata.time,data(i).rawdata(:,2:5));
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
end


plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Frequency (f)')
ylabel('|P1(f)|')

fft_data = fft(data(29).rawdata(:,4));
Fs = 5000;  % sampling frequency
L = 20000;  % length of the signal
f = Fs*(0:(L/2))/L;  % frequency vector
P2 = abs(fft_data/L);  % double-sided spectrum
P1 = P2(1:L/2+1);  % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

% Now, to plot the single-sided amplitude spectrum in log scale
figure;
semilogy(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Frequency (f)')
ylabel('|P1(f)| in dB')
%% FFT

%%%Compute the single-sided amplitude spectrum of the signal.
fft_data = fft(data(29).rawdata(:,4));
Fs = 5000;
L = 20000;
f = Fs*(0:(L-1)/2)/L;
P2 = abs(fft_data/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);
%%In the frequency domain, plot the single-sided spectrum. Because the time sampling of the signal is quite short, the frequency resolution of the Fourier transform is not precise enough to show the peak frequency near 4 Hz.

fft_data = fft(data(1).rawdata(:,2:5),[],1);

plot(data(1).catagoriseddata.time, fft_data)
title ("fft data plottedd")




plot(f(1:500),P1(1:500),"-o") 
title("Single-Sided Spectrum of Original Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")


%% Plot Data - time series 
tiledlayout(8,4)
for i=1:size(data,2);
    nexttile;
    plot(data(i).catagoriseddata.time,data(i).rawdata(:,2:5));
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
end
