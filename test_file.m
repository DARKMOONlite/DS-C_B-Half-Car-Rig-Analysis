%% Setup

addpath(genpath("LAB TEST FOLDER"))
addpath(genpath("functions"))
addpath(genpath("simulink"))
addpath(genpath("data"))

%% Extract Data
clear all
data = Extract_Half_Car_Rig_Data()

%% Plot Data
figure('name','measured data');

tiledlayout(2,4);
for i=1:size(data,2)/4
    nexttile;
    plot(data(i).cdata.time,data(i).rawdata(:,2:5))
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
    
end
legend("m1_translation","m1_roll","m2_translation","m3_translation")

figure('name','converted to dof');
tiledlayout(2,4);
for i=1:size(data,2)/4
    nexttile;
    plot(data(i).cdata.time,data(i).cdata.x1,data(i).cdata.time,data(i).cdata.x2,data(i).cdata.time,data(i).cdata.x3,data(i).cdata.time,data(i).cdata.roll)
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
end
%% FFT

figure("name","fft");
tiledlayout(4,4);
for i=1:length(data)/2
    
    
    nexttile;
    dt=(data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(i).cdata.x1);
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
%% node graphs and tables
figure('name','FFT with nodes');
tiledlayout(4,4)

for i=1:length(data)/2


    nexttile;
    dt = (data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(1).cdata.x1);
    df=1/T;

    fft_data = fft(data(i).rawdof(:,2:5),[],1);
    K = length(fft_data)/2;
    fft_mag = abs(fft_data(1:K,:))/K;
    fft_mag(1,:) = fft_mag(1,:)/2;

    freq = 0:df:(K-1)*df;
    semilogy(freq,fft_mag);
    hold on;
    % Find 4 largest peaks and their corresponding indices for columns 1 to 4
    numPeaks = 4; 
    pks = zeros(numPeaks, 4); % Matrix to store the peak values
    locs = zeros(numPeaks, 4); % Matrix to store the peak indices

   for j = 1:4
        [resolved_peaks, indices] = findpeaks(fft_mag(:, j));
        % Sort peaks in descending order and select the largest ones
        [sortedPeaks, sortIndex] = sort(resolved_peaks, 'descend');
        topPeaks = sortedPeaks(1:min(numPeaks,end));
        topIndices = indices(sortIndex(1:min(numPeaks,end)));
        
        pks(1:length(topPeaks), j) = topPeaks;
        locs(1:length(topIndices), j) = topIndices;
        % Plot peaks on graph
        scatter(freq(topIndices), topPeaks, 25, 'r', '*');
    end

    hold off;
    xlim([0 20]);
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
    for j = 1:4
        node_table = array2table([freq(locs(:,j))' pks(:,j)], 'VariableNames', {'Frequency', 'Magnitude'});
        writetable(node_table,"data/Node_data.xlsx", "Sheet",strcat(data(i).Lift_Position,"_",data(i).Test,"_",data(i).damping,"_",data(i).mass),"Range",strcat(char(j*3+64),string(1)));
    end
end



     %% Create node tables

figure('name','FFT with nodes and phases');
tiledlayout(8,4)  % Double the number of rows to fit phase plots

for i=1:length(data)/2
    nexttile;
    dt = (data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(1).cdata.x1);
    df=1/T;

    fft_data = fft(data(i).rawdof(:,2:5),[],1);
    K = length(fft_data)/2;
    fft_mag = abs(fft_data(1:K,:))/K;
    fft_mag(1,:) = fft_mag(1,:)/2;

    freq = 0:df:(K-1)*df;
    semilogy(freq,fft_mag);
    hold on;
    
   fft_phase = angle(fft_data(1:K,:));  % Phase data
    peakFreqs = zeros(numPeaks, 4);  % Matrix to store the peak frequencies
    peakPhases = zeros(numPeaks, 4);  % Matrix to store the peak phases
       


    for j = 1:4

        [resolved_peaks, indices] = findpeaks(fft_mag(:, j));

        % Sort peaks in descending order and select the largest ones
        [sortedPeaks, sortIndex] = sort(resolved_peaks, 'descend');
        topPeaks = sortedPeaks(1:min(numPeaks,end));
        topIndices = indices(sortIndex(1:min(numPeaks,end)));

        peakFreqs(1:length(topIndices), j) = freq(topIndices);
        peakPhases(1:length(topIndices), j) = fft_phase(topIndices,j);

        freqPhaseTable = table(peakFreqs(:,j), peakPhases(:,j), 'VariableNames', {'Frequency', 'Phase'});
        writetable(freqPhaseTable,"data/Freq_Phase_Data.xlsx", "Sheet",strcat(data(i).Lift_Position,"_",data(i).Test,"_",data(i).damping,"_",data(i).mass),"Range",strcat(char(j*3+64),string(1)));
    end
      % strcat('FreqPhaseTable_', num2str(i), '_column_', num2str(j), '.csv'));
        
        pks(1:length(topPeaks), j) = topPeaks;
        locs(1:length(topIndices), j) = topIndices;
        % Plot peaks on graph
        scatter(freq(topIndices), topPeaks, 25, 'r', '*');
end


hold off;
title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));

% Phase plot
nexttile;
fft_phase = angle(fft_data(1:K,:));
plot(freq, rad2deg(fft_phase));  % Convert to degrees
hold on;
for j = 1:4
    scatter(freq(locs(:,j)), rad2deg(fft_phase(locs(:,j),j)), 25, 'r', '*');
end
hold off;
title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass," Phase"));


%% Waterfall plot 

% Assuming freq and fft_mag are your frequency and FFT magnitude data respectively.
% freq should be a 1D array, while fft_mag should be a 2D array, with each column representing data from each run
resolution = 1;
i=1;
ribbon_data = zeros(length(fft(data(i).rawdof(:,2:5),[],1))/2,4)
for i =1:4
    dt = (data(i).cdata.time(2)-data(i).cdata.time(1));
    T = dt*length(data(1).cdata.x1);
    df=1/T;

    fft_data = fft(data(i).rawdof(:,2:5),[],1);
    K = length(fft_data)/2;
    fft_mag_w = abs(fft_data(1:K,:))/K;
    fft_mag_w(1,:) = fft_mag_w(1,:)/2;
    ribbon_data(:,i) = fft_mag_w(1:resolution:end,1);

end


waterfall_temp = zeros(size(ribbon_data));
for i=1:length(waterfall_temp)
    waterfall_temp(i,:)=[500,1000,1500,2000];
end
waterfall_time=zeros(size(ribbon_data));
for i=1:size(waterfall_temp,2)
    waterfall_time(:,i) = [0:df:df*(length(waterfall_temp)-1)];
end

waterfall(waterfall_temp(1:K/50,:),waterfall_time(1:K/50,:),mag2db(ribbon_data(1:K/50,:)));

% semilogy(waterfall_time(:,1),ribbon_data(:,1))
% axis([0 100 1e-6 10]);

%% TRANSFER FUNCTION - Welch and Hann comparison to fft, undamped no mass 
init_var = data(29).init_var;
fft_data = fft(data(29).rawdof(:,2), [], 1);

% Ensure init_var and fft_data have the same size
if size(init_var,1) < size(fft_data,1)
    init_var = repmat(init_var,size(fft_data,1),1);
elseif size(init_var,1) > size(fft_data,1)
    fft_data = repmat(fft_data,size(init_var,1),1);
end

Pack1 = [init_var fft_data];

blocks = 1;
NFFT = 20000 * blocks;
fs = 5000; % Define sampling frequency
count = 1;

% Here you need to define which column from Pack1 is the input and which is the output
input = Pack1(:,1);
output = Pack1(:,2);

[Gd, F] = cpsd(input, input, hann(NFFT), NFFT/2, NFFT, fs);
[Gn, ~] = cpsd(output, input, hann(NFFT), NFFT/2, NFFT, fs);

%[Gd, F] = cpsd(input, input, hann(NFFT), NFFT/2, NFFT, fs); % PSD of input
%[Gn, ~] = cpsd(output, output, hann(NFFT), NFFT/2, NFFT, fs); % PSD of output
%[Gcross, ~] = cpsd(input, output, hann(NFFT), NFFT/2, NFFT, fs); % CPSD of input and output
%H(:, count) = Gcross ./ Gd; % Transfer function is the ratio of the CPSD and the PSD of input

H(:, count) = Gn ./ Gd;
count = count + 1;
  
H1A = mean(H, 2);
H1A_mag = abs(H1A);
H1A_mag = H1A_mag(1:length(F)); % Match the size with F
H1A_phase = -angle(H1A) * 180/pi;
H1A_phase = H1A_phase(1:length(F)); % Match the size with F

clear H Gn Gd count;

% FRF
figure;
subplot(2, 1, 1);
    plot(F, 20*log10(H1A_mag));
    grid on; grid minor;
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title('Center undamped no mass lab data transfer function using Welch and Hann window');
    xlim([0 10]);

subplot(2, 1, 2);
    plot(F, H1A_phase);
    grid on; grid minor;
    xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
    xlim([0 10]);

    %%phase not working input = Pack1(:,1);

    %% Now fit the transfer function 
   dataWTF = iddata(output, input, 1/1000);
    symWTF = tfest(dataWTF,3)
    sysdWTF = c2d(symWTF,1/1000);
    estimate1 = filter(sysdWTF.Numerator,sysdWTF.Denominator, init_var);

