%% Extract Data
data = Extract_Half_Car_Rig_Data()

%% Plot Data
tiledlayout(8,4)
for i=1:size(data,2)
    nexttile
    plot(data(i).catagoriseddata.time,data(i).rawdata(:,2:5))
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass))
end

%% FFT


fft_data = fft(data(1).rawdata(:,2:5),[],1)

plot(data(1).catagoriseddata.time, fft_data);

%%%Compute the single-sided amplitude spectrum of the signal.
fft_data = fft(data(29).rawdata(:,4));
Fs = 5000;
L = 20000;
f = Fs*(0:(L-1)/2)/L;
P2 = abs(fft_data/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);
%%In the frequency domain, plot the single-sided spectrum. Because the time sampling of the signal is quite short, the frequency resolution of the Fourier transform is not precise enough to show the peak frequency near 4 Hz.

%%plot phase angle 
phase = -angle(P2)*180/pi;
plot (f(1:500),phase(1:500),"-o")
grid on;grid minor
xlabel("f (Hz)")
ylabel("phase")


plot(f(1:500),P1(1:500),"-o") 
title("Single-Sided Spectrum of Original Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")
