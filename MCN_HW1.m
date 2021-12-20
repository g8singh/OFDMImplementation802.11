close all
clear 
clc

%-----------------------------------------------------------------------
% (1) Packet construction and OFDM modulation.
%-----------------------------------------------------------------------

% Generate data of 4160 bits by random manner

org_data = randi([0 1],4160,1);
data = [reshape(org_data,1,4160) zeros(1,16)];

% Use the Pilot Data

pilot = zeros(4,1);

% Use the BPSK Modulation and Demodulation at transmitter and receiver
% respectively

bModulator   = comm.BPSKModulator;
bDeModulator = comm.BPSKDemodulator;

N       = 64;

N_SD    = 48;
N_SP    = 4;
N_ST    = N_SD + N_SP;

Delta_F = 0.3125*(10^6);

format long
T_FFT = 1/Delta_F;

T_Short = 10*(T_FFT/4);
T_GI    = (T_FFT/4);
T_GI2   = 2*T_GI;
T_Long  = (T_GI2 + (2*T_FFT));
T_Sym   = (T_GI + T_FFT);

N_CP    = ((N*T_GI)/T_FFT);

% Generate the BPSK Signal

modData      = bModulator(reshape(data,4176,1));
pilotData    = bModulator(pilot);

modData      = reshape(modData, 1, 4176);

stfData      = [0 0 complex(1,1) 0 0 0 complex(-1,-1) 0 0 0 complex(1,1) 0 0 0 complex(-1,-1) 0 0 0 complex(-1,-1) 0 0 0 complex(1,1) 0 0 0 0 0 0 0 complex(-1,-1) 0 0 0 complex(-1,-1) 0 0 0 complex(1,1) 0 0 0 complex(1,1) 0 0 0 complex(1,1) 0 0 0 complex(1,1) 0 0];
stfData      = sqrt(13/6)*stfData;

ltfData      = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];

% Use OFDM Demodulation

X_n          = zeros(87,64);
X_t          = zeros(87,64);
X_t_ofdm     = zeros(87,80);

for i = 1:87
    k_offset     = (i-1)*48;
    k_lower      = k_offset + 1;
    k_upper      = k_lower + 24;
    format long
    %X_n(i,:)     = [zeros(1,6) pilotData(1) modData(k_lower:k_lower+23) pilotData(2) zeros(1,1) pilotData(3) modData(k_upper:k_upper+23) pilotData(4) zeros(1,5)];
    X_n(i,:)     = [zeros(1,1) pilotData(1) modData(k_lower:k_lower+23) pilotData(2) zeros(1,11) pilotData(3) modData(k_upper:k_upper+23) pilotData(4)];
    %X_n(i,:)     = [zeros(1,1) modData(k_lower_1:k_lower_1+5) pilotData(1) modData(k_lower_2:k_lower_2+12) pilotData(2) modData(k_lower_3:k_lower_3+4) zeros(1,11) modData(k_upper_1:k_upper_1+4) pilotData(3) modData(k_upper_2:k_upper_2+12) pilotData(4) modData(k_upper_3:k_upper_3+5)];
    format long
    X_t(i,:)     = ifft(X_n(i,:), N);
    X_t_cp       = X_t(i, 49:64);
    X_t_ofdm(i,:)= [X_t_cp X_t(i,:)];
end

% Generate packet preambles

%Stf_n        = [zeros(1,6) stfData(1:53) zeros(1,5)];
format long
Stf_n        = [stfData(27) stfData(28:53) zeros(1,11) stfData(1:26)];

format long
Ltf_n        = [ltfData(27) ltfData(28:53) zeros(1,11) ltfData(1:26)];

format long
stf_t        = ifft(Stf_n, 64);
format long
ltf_t        = ifft(Ltf_n, 64);

%Plot the PSD and the STF magnitude

for i = 1:height(X_t_ofdm)
    xdft         = fft(X_t_ofdm(i,:),N);
    psd          = (1/(N)) * abs(xdft).^2;
    freq         = 1:1:N;
    if (i == 1)
        figure
    elseif(i == 2)
        hold on
    end
    plot(freq,psd)
end
title('Power Spectral Density of OFDM Symbols Using FFT')
xlabel('Frequency') 
ylabel('Power/Frequency')
hold off
savefig('PSD.fig')
close(figure)

stf_t_ofdm = [stf_t stf_t stf_t(1:32)];

figure
stem(abs(stf_t_ofdm))
title('Plot of STF magnitude before the channel distortion')
xlabel('Sample') 
ylabel('Magnitude')
savefig('STF_Mag_Before_Channel.fig')
close(figure)

ltf_t_ofdm = [ltf_t(33:64) ltf_t ltf_t];

x_t_ofdm = [];

% Serialize the packets into a single dimensional array

for i = 1:height(X_t_ofdm)
    x_t_ofdm = [x_t_ofdm X_t_ofdm(i,:)];
end

%-----------------------------------------------------------------------
% (2) Packet transmission and channel distortion
%-----------------------------------------------------------------------

% Add the idle sequence
fin_signal = [zeros(1,100) stf_t_ofdm ltf_t_ofdm x_t_ofdm];

% Add the channel attenuation and phase shift
fin_signal_Amp_Mod = fin_signal*(10^(-5));
fin_signal_Phase_Shift = fin_signal_Amp_Mod*exp(complex(0,((-3*pi)/4)));

fin_signal_freq_offset_noise = zeros(1,length(fin_signal_Phase_Shift));

% Add the frequency offset and noise
for i = 1:length(fin_signal_Phase_Shift)
    k = (i-1);
    fin_signal_freq_offset_noise(i) = fin_signal_Phase_Shift(i)*exp(complex(0,(-2*pi*0.00017*k)));
    %fin_signal_freq_offset_noise(i) = fin_signal_freq_offset_noise(i) + complex((randn(1)*(10^(-11))),(randn(1)*(10^(-11))));
end

fin_signal_freq_offset_noise = awgn(fin_signal_freq_offset_noise,20,'measured');

%Plot the STF magnitude after channel distortion
figure
stem(abs(fin_signal_freq_offset_noise(101:260)))
title('Plot of STF magnitude after the channel distortion')
xlabel('Sample') 
ylabel('Magnitude')
savefig('STF_Mag_After_Channel.fig')
close(figure)

%-----------------------------------------------------------------------
% (3) Packet detection
%-----------------------------------------------------------------------

% Generating the self-correlation function(R(m)) and the energy
% function(E(m))
R_t = zeros(1,length(fin_signal_freq_offset_noise));
E_t = zeros(1,length(fin_signal_freq_offset_noise));

first_idx = 0;

for m = 17:(length(fin_signal_freq_offset_noise)-16)
    sum_R = 0;
    sum_E = 0;
    for i = 0:15
        k_ref = (m+i);
        sum_R = sum_R + fin_signal_freq_offset_noise(k_ref)*conj(fin_signal_freq_offset_noise(k_ref-16));
        sum_E = sum_E + fin_signal_freq_offset_noise(k_ref)*conj(fin_signal_freq_offset_noise(k_ref));
    end
    R_t(m-16) = sum_R;
    E_t(m-16) = sum_E;
    if((sum_R >= 5*(10^(-13))) && (first_idx == 0))
        first_idx = m;
    end
end

% Plot the self-correlation function and energy function
figure
t = 1:1:length(R_t);
plot(t, (abs(R_t)))
hold on
plot(t, (abs(E_t)))
hold off
legend('R(m)','E(m)')
title('Plot of self correlation values(R(m)) and signal energy(E(m))')
xlabel('Sample(m)') 
ylabel('Magnitude of the values')
savefig('Self_Correlation.fig')
close(figure)

%-----------------------------------------------------------------------
% (4) Packet synchronization
%-----------------------------------------------------------------------

% Generating the cross-correlation function R_t(m)
for m = 1:(length(fin_signal_freq_offset_noise)-16)
    sum_R = 0;
    for i = 0:15
        k_ref = (m+i);
        sum_R = sum_R + fin_signal_freq_offset_noise(k_ref)*conj(stf_t((i+1)));
    end
    R_t(m) = sum_R;
end

% Plotting the cross-correlation function
figure
plot(t, abs(R_t))
title('Plot of correlation values(R(m)) with STF signal')
xlabel('Sample(m)') 
ylabel('Magnitude of the values')
savefig('STF_Correlation.fig')
close(figure)

% Pick the indexes with 10 peaks
[pks, locs] = maxk(abs(R_t), 10);

locs = sort(locs);

%disp(locs)

%-----------------------------------------------------------------------
% (5) Channel Estimation and Packet decoding
%-----------------------------------------------------------------------

%First instance of the peak is where STF starts. 
first_stf_indx    = locs(1);
first_ltf_indx    = first_stf_indx + length(stf_t_ofdm);
first_ltf_indx_s1 = first_ltf_indx + 32;

delta_f = zeros(1,(length(ltf_t_ofdm)-N));

%Estimate the frequency offset through LTF
for i = first_ltf_indx:(first_ltf_indx + length(ltf_t_ofdm) - N - 1)
    delta_f((i-first_ltf_indx)+1) = (imag(fin_signal_freq_offset_noise(i+N)/fin_signal_freq_offset_noise(i)))/(2*pi*N);
end

delta_f_avg = sum(delta_f)/length(delta_f);

disp(delta_f_avg)

dec_signal_after_freq_offset = zeros(1,length(fin_signal_freq_offset_noise));

%Add the estimated frequency offset to the received signal
for i = 1:length(fin_signal_freq_offset_noise)
    k = (i-1);
    dec_signal_after_freq_offset(i) = fin_signal_freq_offset_noise(i)*exp(complex(0,(-2*pi*delta_f_avg*k)));
    %fin_signal_freq_offset_noise(i) = fin_signal_freq_offset_noise(i) + (randn(1)*(10^(-11)));
end

%Determine the channel gain
Y_k = fft(dec_signal_after_freq_offset(first_ltf_indx_s1:first_ltf_indx_s1+N-1),N);

H_k = zeros(1,N);

for i = 1:N
    if(Ltf_n(i) ~= 0)
        H_k(i) = Y_k(i)/Ltf_n(i);
    else
        H_k(i) = 0;    
    end
end

disp(abs(H_k))
disp(angle(H_k)/pi)

figure
stem(abs(H_k))
title('Plot of magnitude of channel again wrt frequency')
ylabel('Magnitude') 
xlabel('Frequency')
savefig('Hk_Mag.fig')
close(figure)

figure
stem(angle(H_k)/pi)
title('Plot of phase of channel again wrt frequency')
ylabel('Phase(In multiples of \pi)') 
xlabel('Frequency')
savefig('Hk_Phase.fig')
close(figure)


%Get the start index of the data symbols by using the known LTF length and
%start index
ofdm_data_start_idx = first_ltf_indx + length(ltf_t_ofdm);

X_decoded = zeros(1,N);
X_conc    = [];

% Generate the expected sequence after performing FFT and dividing with channel gain
for i = ofdm_data_start_idx:80:(length(dec_signal_after_freq_offset)-80+1)
    X_k = fft(dec_signal_after_freq_offset((i+16):(i+79)), N);
    for j = 1:N
        if(H_k(j) ~= 0)
            %disp(angle(X_k(j)) - angle(H_k(j)))
            X_decoded(j) = X_k(j)/H_k(j);
        else
            X_decoded(j) = 0;
        end
    end
    X_conc = [X_conc X_decoded(3:26) X_decoded(40:63)];
    %X_conc = [X_conc X_decoded(2:7) X_decoded(9:21) X_decoded(23:26) X_decoded(39:43) X_decoded(45:57) X_decoded(59:64)];
end

X_dec = zeros(1,(length(X_conc)-16));

%disp(angle(X_conc))

% Use phase approximation to get the right BPSK symbols
for i = 1:length(X_conc)-16
    if((angle(X_conc(i)) > (pi/2)) || (angle(X_conc(i)) < -(pi/2)))
        X_dec(i) = -1;
    else
        X_dec(i) = 1;
    end
end

% Demodulate the BPSK modulated signal
X_bin = bDeModulator(reshape(X_dec,length(X_dec),1));

num_errors = 0;

error_plot = zeros(1,length(X_bin));

for i = 1:length(X_bin)
    if(X_bin(i) ~= org_data(i))
        num_errors = num_errors + 1;
    end
    disp([X_bin(i) org_data(i)])
    error_plot(i) = num_errors;
end

figure
stem(error_plot)
title('Plot of cumulative error wrt data index')
xlabel('Data index') 
ylabel('Cumulative error')
savefig('Error_Plot.fig')
close(figure)
