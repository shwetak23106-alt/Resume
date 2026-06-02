%==========================================================================
% TITLE:  Biomedical Signal Processing: 50 Hz Powerline Interference Removal
% AUTHOR: Shweta K
% DATE:   June 2026
% DATA SOURCE: PhysioNet MIT-BIH Arrhythmia Database (Record 100)
% DESCRIPTION: This script unpacks raw, compressed Format 212 ECG data
%              and applies a high-selectivity, recursive 2nd-order IIR 
%              Notch Filter with zero-phase correction to suppress 
%              environmental powerline hum while preserving P-Q-R-S-T morphology.
%==========================================================================

%% --- STAGE 1: RAW BINARY EXTRACTION & CALIBRATION (FORMAT 212 DECODER) ---
clear; clc; close all;

Fs = 360;           % Native sampling rate of the MIT-BIH database (360 samples/sec)
num_samples = 3000; % Time-window boundary: Extract first ~8.3 seconds of continuous data

% Establish a low-level file identifier link with the binary dataset
fid = fopen('100.dat', 'r');
if fid == -1
    error('IO_ERROR: 100.dat not found. Ensure raw PhysioNet files are in the active directory.');
end

% MITDB Format 212 Packing: Read raw data stream sequentially as unsigned 8-bit blocks (bytes)
raw_bytes = fread(fid, [3, num_samples], 'uint8')';
fclose(fid);

% Bitwise manipulation array initialization for Channel 1 (MLII Lead)
ecg_channel1 = zeros(num_samples, 1);

for i = 1:num_samples
    % DEMULTIPLEXING LOGIC: Isolate the lower 8 bits from Byte 1 and combine them with 
    % the lower 4-bit nibble of Byte 2 (isolated via 15/0x0F mask) shifted left to the MSB position.
    ecg_channel1(i) = bitand(raw_bytes(i,1), 255) + bitshift(bitand(raw_bytes(i,2), 15), 8);
end

% CALIBRATION: Convert unscaled digital ADC values to physical voltage magnitudes (mV)
% by subtracting the standard 1024 digital baseline offset and dividing by the 200x amplifier gain factor.
real_patient_ecg = (ecg_channel1 - 1024) / 200; 
t = (0:num_samples-1) / Fs; % Map the continuous temporal axis vector based on Fs

%% --- STAGE 2: MATHEMATICAL ARTIFACT INJECTION ---
% Simulate localized electromagnetic powerline coupling by superimposing a heavy 
% 50 Hz pure sinusoidal wave (0.4 mV amplitude) onto the clean baseline clinical signal.
noise_50Hz = 0.4 * sin(2 * pi * 50 * t)'; 
noisy_ecg = real_patient_ecg + noise_50Hz;

%% --- STAGE 3: DIGITAL FILTER DESIGN & ZERO-PHASE PROCESSING ---
% Normalize target cutoff frequency relative to the discrete system's Nyquist boundary (Fs/2 = 180 Hz)
wo = 50 / (Fs/2);                     

% Define notch selectivity: Q-Factor approximation establishing highly narrow stopband bandwidth
bw = wo / 30;                         

% Compute 2nd-order feedforward (b) and feedback (a) transfer function coefficients via Z-plane mapping
[b, a] = iirnotch(wo, bw);            

% Apply a 2nd-order recursive IIR Notch Filter via pole-zero placement to surgically 
% attenuate the 50Hz powerline hum while keeping computational overhead minimal.
% Dual-pass execution (forward and backward) neutralizes phase delays to preserve peak timing.
cleaned_ecg = filtfilt(b, a, noisy_ecg); 

%% --- STAGE 4: SPECTRUM ANALYSIS (DISCRETE FOURIER TRANSFORM) ---
N = length(t);
f = (0:N-1) * (Fs / N);               % Frequency scaling vector spanning the discrete spectrum

% Compute Fast Fourier Transform (FFT) magnitudes normalized by vector length N 
% to extract real-world spectral voltage densities
FFT_noisy = abs(fft(noisy_ecg)) / N;
FFT_cleaned = abs(fft(cleaned_ecg)) / N;

%% --- STAGE 5: VISUALIZATION DASHBOARD GENERATION ---
figure('Name', 'PhysioNet ECG Processing Dashboard', 'NumberTitle', 'off');

% Plot 1: Time-Domain Corrupted Input
subplot(2,2,1); plot(t, noisy_ecg, 'r');
title('Real Patient ECG + Injected 50Hz Noise');
xlabel('Time (seconds)'); ylabel('Amplitude (mV)'); grid on;

% Plot 2: Time-Domain Recovered Waveform
subplot(2,2,3); plot(t, cleaned_ecg, 'g');
title('Filtered Output (IIR Notch Cleaned)');
xlabel('Time (seconds)'); ylabel('Amplitude (mV)'); grid on;

% Plot 3: Frequency-Domain Input Profile
subplot(2,2,2); plot(f, FFT_noisy, 'r');
xlim([0 100]); title('Frequency Spectrum (Before Filter)');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;

% Plot 4: Frequency-Domain Output Profile
subplot(2,2,4); plot(f, FFT_cleaned, 'g');
xlim([0 100]); title('Frequency Spectrum (After Filter)');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;