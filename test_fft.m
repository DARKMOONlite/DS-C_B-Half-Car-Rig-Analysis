% clear all; close all; clc
% %% load data
% x200 = load('200_rect_cut.csv');
% %x200H = load('200_Hann_cut.csv'); %not required since time data are the
% %same in each dataset
% x199 = load('199_rect_cut.csv');
%x199H = load('199_Hann_cut.csv'); %ditto as above
%% extract sampling parameters
dt = x200(2,1)-x200(1,1); %time between samples
Fs = 1/dt; %sample frequency
T = dt*length(x200(:,2)); %determine the no. of samples in the time array
df = 1/T;
%% process data
%create ffts
%200
fftx200 = fft(x200(:,2));
K = length(fftx200)/2; %determine size of spectrum (= half length of time array)
xf200mag = sqrt(fftx200(1:K).*conj(fftx200(1:K))); %mutliply up to Nyquist bycomplex conj. to get magnitude spectrum
xf200mag = xf200mag*2; %double the amplitude to account for folding spectrum
xf200mag(1) = xf200mag(1)/2; %fix the DC component which should not be doubled
xf200mag = xf200mag/length(x200); xf200mag = xf200mag'; %fix the scaling to account for the length
%apply Hann
fftx200H = fft(x200(:,2).*hann(length(x200(:,2)))); %multiply input time array by Hann window
xf200Hmag = sqrt(fftx200H(1:K).*conj(fftx200H(1:K)));
xf200Hmag = xf200Hmag*2;
xf200Hmag(1) = xf200Hmag(1)/2;
xf200Hmag = xf200Hmag/length(x200); xf200Hmag=xf200Hmag';
%199
fftx199 = fft(x199(:,2));
xf199mag = sqrt(fftx199(1:K).*conj(fftx199(1:K)));
xf199mag = xf199mag*2;
xf199mag(1) = xf199mag(1)/2;
xf199mag = xf199mag/length(x199); xf199mag=xf199mag';
%apply Hann
fftx199H = fft(x199(:,2).*hann(length(x199(:,2))));
xf199Hmag = sqrt(fftx199H(1:K).*conj(fftx199H(1:K)));
xf199Hmag = xf199Hmag*2;
xf199Hmag(1) = xf199Hmag(1)/2;
xf199Hmag = xf199Hmag/length(x199); xf199Hmag=xf199Hmag';
%% plot the spectra
%create freq vector for the x-axis
freq = 0:df:(K-1)*df; %only need to do this once for al spectra
figure(1)
subplot(2,1,1)
semilogy(freq,xf200mag); axis([0 500 1e-5 10]); title('200 Hz rectangular window')
hold on
subplot(2,1,2)
semilogy(freq,xf200Hmag); axis([0 500 1e-5 10]); title('200 Hz Hann window')
figure(2)
subplot(2,1,1)
semilogy(freq,xf199mag); axis([0 500 1e-5 10]); title('199 Hz rectangular window')
hold on
subplot(2,1,2)
semilogy(freq,xf199Hmag); axis([0 500 1e-5 10]); title('199 Hz Hann window')
