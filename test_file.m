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



%%ribbon
filtered_data = data(indices);
% Initialize an empty array to store indices
indices = [];

% Loop through the data structure
for i = 1:length(data)
    % Check for conditions in the desired columns
    if strcmp(data(i).rawdof(4,:),'damped') && strcmp(data(i).rawdof(6,:),'centre')
        % If conditions are met, append the index to the array
        indices = [indices; i];
    end
end

% Use indices to get the desired data from the structure
fft_results = cell(size(filtered_data));  % initialize an empty cell array to store fft results

% Apply FFT
for i = 1:length(filtered_data)
    fft_results{i} = fft(filtered_data(i).rawdof(:,2),[],1);
end
% Convert the cell array to a matrix
fft_matrix = cell2mat(fft_results');

% Create the ribbon plot
figure;
ribbon(abs(fft_matrix)); % abs is used to get the magnitude of the complex numbers resulted from FFT
title('Ribbon Plot');

%%%
fft(data(i).rawdof(:,2:5),[],1);
ribbon = fft(data(i).rawdof(:,2),[],1)


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



% Make sure freq is a column vector
% freq = freq(:);
% 
% % Make sure fft_mag is oriented correctly
% if size(fft_mag, 1) ~= length(freq)
%     fft_mag = fft_mag.'
% end
% 
% % Create a grid of y values (each run)
% Y = 1:size(fft_mag, 2)
% 
% % Create a grid for the waterfall plot
% [X, Y] = meshgrid(freq, Y)
% 
% % Plot the waterfall
% figure
% waterfall(X, Y, fft_mag.')
% xlabel('Frequency (Hz)')
% ylabel('Run Number')
% zlabel('FFT Magnitude')
% title('Waterfall Plot')
% 
% 
% ribbon
% % meshgrid(-3:.125:3);
% 
% %%gm FFT 
% % Extract Data
% %data = Extract_Half_Car_Rig_Data();
% 
% Fs = 5000;  % sampling frequency
% 
% % Create tiled layout
% tiledlayout(8,4)
% 
% % Iterate over all data sets
% for i=1:size(data,2)
%     % Compute FFT of the data
%     fft_data = fft(data(i).rawdata(:,4));
%     L = length(data(i).rawdata(:,4));  % length of the signal
%     f = Fs*(0:(L/2))/L;  % frequency vector
%     P2 = abs(fft_data/L);  % double-sided spectrum
%     P1 = P2(1:L/2+1);  % single-sided spectrum
%     P1(2:end-1) = 2*P1(2:end-1);
% 
%     % Plot the single-sided amplitude spectrum in the next tile
%     nexttile;
%     semilogy(f,P1)
%     title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass));
%     xlabel('Frequency (f)')
%     ylabel('|P1(f)| in dB')
% 
%     % Set the x-axis limit to 20
%     xlim([0 20])
% end


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

%% % find peaks and their corresponding indices
    
    [pks, locs] = findpeaks(fft_mag (:,1), 'MinPeakHeight', 0.0001, 'MinPeakDistance', 0.001);
    