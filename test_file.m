

%% Extract Data
clear all
data = Extract_Half_Car_Rig_Data()

%% Plot Data
figure('name','measured data')

tiledlayout(8,4);
for i=1:size(data,2)
    nexttile;
    plot(data(i).cdata.time,data(i).rawdata(:,2:5))
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass))
end

figure('name','converted to dof');
tiledlayout(8,4);
for i=1:size(data,2)
    nexttile;
    plot(data(i).cdata.time,data(i).cdata.x1,data(i).cdata.time,data(i).cdata.x2,data(i).cdata.time,data(i).cdata.x3,data(i).cdata.time,data(i).cdata.roll)
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass))
end
%% FFT

figure("name","fft");
tiledlayout(4,4);
for i=1:length(data)/2
    
    
    nexttile;
    dt=(data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(1).cdata.x1);
    df=1/T;

    fft_data= fft(data(i).rawdof(:,2:5),2^nextpow2(length(data(i).cdata.time)));
    K = length(fft_data)/2;
    fft_mag = sqrt(fft_data(1:K,:).*conj(fft_data(1:K,:)));
    fft_mag = fft_mag*2;
    fft_mag(1)=fft_mag(1)/2;
    fft_mag = fft_mag/length(fft_mag);
    
    % bins = [0:K-1]
    freq = 0:df:(K-1)*df;

    semilogy(freq,fft_mag);
    axis([0 100 1e-6 10]);
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
end
%% these ffts have all been multiplied so that their highest amplitued in time domain is 1
figure("name"," multiplied fft");
tiledlayout(4,4);
for i=1:length(data)/2
    init_var = data(i).init_var;
    multiplied_data = data(i).rawdof(:,:);
    for j=1:length(init_var)
        multiplier(j) = 1000/init_var(j);
        multiplied_data(:,j+1) =  multiplied_data(:,j+1)*multiplier(j);
    end
    

    
    

    nexttile;
    dt=(data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(1).cdata.x1);
    df=1/T;

    fft_data= fft(multiplied_data(:,2:5),2^nextpow2(length(data(i).cdata.time)));
    K = length(fft_data)/2;
    fft_mag = sqrt(fft_data(1:K,:).*conj(fft_data(1:K,:)));
    fft_mag = fft_mag*2;
    fft_mag(1)=fft_mag(1)/2;
    fft_mag = fft_mag/length(fft_mag);
    
    % bins = [0:K-1]
    freq = 0:df:(K-1)*df;

    semilogy(freq,fft_mag);
    axis([0 100 1e-6 10]);
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
end


% fft_data = fft(data(1).rawdata(:,2:5),[],1)
% 
% plot(data(1).cdata.time, fft_data)

%%%Compute the single-sided amplitude spectrum of the signal.
% fft_data = fft(data(1).rawdata(:,2:5),[],1);
% Fs = 5000;
% L = 20000;
% f = Fs*(0:(L-1)/2)/L;

% dt=(data(1).cdata.time(2)-data(1).cdata.time(1));
% T = dt*length(data(1).cdata.LVDT1);
% df=1/T;
% 
% K = length(fft_data)/2;
% 
% fft_mag = sqrt(fft_data(1:K,1).*conj(fft_data(1:K,1)));
% fft_mag = fft_mag*2;
% fft_mag(1)=fft_mag(1)/2;
% fft_mag = fft_mag/length(fft_mag);
% 
% 
% 
% freq = 0:df:(K-1)*df
% 
% 
% semilogy(freq,fft_mag)
% axis([0 500 1e-5 10]); title('200 Hz rectangular window');
% P2 = abs(fft_data/L);
% P1 = P2(1:(L+1)/2);
% P1(2:end) = 2*P1(2:end);
% %%In the frequency domain, plot the single-sided spectrum. Because the time sampling of the signal is quite short, the frequency resolution of the Fourier transform is not precise enough to show the peak frequency near 4 Hz.
% 
% plot(f,P1,"-o") 
% title("Single-Sided Spectrum of Original Signal")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")
