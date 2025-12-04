Fs = 8000;


recObj = audiorecorder(Fs, 16, 1);
disp('Ξεκινά η ηχογράφηση... Πείτε το μικρό σας όνομα (διάρκεια < 3 sec).');
recordblocking(recObj, 2.5);
disp('Τέλος ηχογράφησης.');
x = getaudiodata(recObj);
m_path = fileparts(mfilename('fullpath'));
filename = 'panagiotis.wav'; 
output_path = fullfile(m_path, filename);
audiowrite(output_path, x, Fs);



[x, Fs] = audioread('panagiotis.wav'); 


t = (0:length(x)-1)/Fs;

% Σχεδίαση στο πεδίο του χρόνου
figure('Name', 'Ερώτημα 1α - Σήμα Φωνής');
plot(t, x);
title('Σήμα Φωνής στο Πεδίο του Χρόνου');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος');
grid on;


x_norm = x / max(abs(x));


M = 400; 
x_sq = x_norm.^2; 


w = hamming(M);


Energy = conv(x_sq, w);


%% 1γ) Υπολογισμός Zero Crossing Rate (ZCR)



s = sign(x_norm);

zero_crossings_signal = abs(s(2:end) - s(1:end-1));

zero_crossings_signal(zero_crossings_signal > 1) = 1; 
zero_crossings_signal = [0; zero_crossings_signal];


w_zcr = (1/(2*M)) * ones(M, 1);


ZCR = conv(zero_crossings_signal, w_zcr);


%% 1δ) Σχεδίαση σε Κοινό Διάγραμμα και Σχολιασμός


L = length(x_norm);



Energy_scaled = Energy(floor(M/2):L + floor(M/2) - 1);
ZCR_scaled = ZCR(floor(M/2):L + floor(M/2) - 1);


Energy_scaled = Energy_scaled / max(Energy_scaled);
ZCR_scaled = ZCR_scaled / max(ZCR_scaled);


figure('Name', 'Ερώτημα 1δ - Σήμα, Ενέργεια και ZCR');
plot(t, x_norm, 'b', 'DisplayName', 'Σήμα Φωνής x[n]');
hold on;
plot(t, Energy_scaled, 'r', 'DisplayName', 'Ενέργεια');
plot(t, ZCR_scaled, 'g', 'DisplayName', 'ZCR');
title('Σήμα Φωνής, Ενέργεια και ZCR');
xlabel('Χρόνος (sec)');
ylabel('Κλιμακωμένο Πλάτος');
legend('Location', 'best');
grid on;

% [ΣΧΟΛΙΑΣΜΟΣ]
% (Γράψτε στην αναφορά σας: Παρατηρείται ότι οι έμφωνοι ήχοι (υψηλό πλάτος/ενέργεια) 
% έχουν χαμηλό ZCR, ενώ οι άφωνοι ήχοι (χαμηλό πλάτος/ενέργεια) έχουν υψηλό ZCR. 
% Οι μετρήσεις αυτές είναι κρίσιμες για το διαχωρισμό τους.)


%% 1ε) Επιλογή Περιοδικού Αποσπάσματος και Υπολογισμός Περιόδου



L_seg = 0.05 * Fs; 


n_start = 1000; 
x_seg = x_norm(n_start : n_start + L_seg - 1);
t_seg = (0:L_seg-1)/Fs;


disp('Ακρόαση αποσπάσματος 50 msec...');
sound(x_seg, Fs);


figure('Name', 'Ερώτημα 1ε - Απόσπασμα 50 msec');
plot(t_seg, x_seg);
title(['Απόσπασμα 50 msec (', num2str(L_seg), ' δείγματα)']);
xlabel('Χρόνος (sec)');
ylabel('Πλάτος');
grid on;

% [ΣΧΟΛΙΑΣΜΟΣ]
% (Γράψτε στην αναφορά σας: Το απόσπασμα αντιστοιχεί στο φώνημα [...]. 
% Μετρώντας την απόσταση μεταξύ δύο διαδοχικών κορυφών βρέθηκε η περίοδος T0 = [...].)


%% 1στ) Διακριτός Μετασχηματισμός Fourier (DFT)


N = 1024; 
Xk = fft(x_seg, N);


Magnitude_Xk = abs(Xk);
Magnitude_dB = 20*log10(Magnitude_Xk);


N_half = N/2;
f = (0:N_half-1) * (Fs/N); 


Magnitude_linear_half = Magnitude_Xk(1:N_half);
Magnitude_dB_half = Magnitude_dB(1:N_half);


figure('Name', 'Ερώτημα 1στ - Φάσμα');


subplot(2,1,1);
plot(f, Magnitude_linear_half);
title('Μέτρο Φάσματος (|X[k]|) - Γραμμική Κλίμακα');
xlabel('Συχνότητα (Hz)');
ylabel('Πλάτος');
grid on;


subplot(2,1,2);
plot(f, Magnitude_dB_half);
title('Μέτρο Φάσματος (20log_{10}|X[k]|) - Λογαριθμική Κλίμακα');
xlabel('Συχνότητα (Hz)');
ylabel('Πλάτος (dB)');
grid on;


%% 1ζ) Υπολογισμός Θεμελιώδους Συχνότητας

[~, peak_index] = max(Magnitude_linear_half(2:end));
F0_index = peak_index + 1; 

F0_Hz = f(F0_index);

disp(['Η Θεμελιώδης Συχνότητα F0 (εποπτικά από το φάσμα) είναι: ', num2str(F0_Hz), ' Hz']);

% [ΣΧΟΛΙΑΣΜΟΣ]
% (Γράψτε στην αναφορά σας: Η εποπτική περίοδος T0 (από 1ε) ήταν [...] sec. 
% Η αντίστοιχη θεμελιώδης συχνότητα F0 = 1/T0 είναι [...] Hz. 
% Αυτή η τιμή επιβεβαιώνεται από την F0 που βρέθηκε στο φάσμα (1ζ) η οποία είναι [...] Hz.)



% ΕΝΟΤΗΤΑ 2: ΣΧΕΔΙΑΣΗ ΦΙΛΤΡΩΝ

%% 2.1α) Σχεδίαση Φίλτρου, Πόλοι-Μηδενικά και Συντελεστές
z = [1.2; -0.6];
p = [0.196 + 0.672i; 0.196 - 0.672i];
K = 0.15;

figure('Name', 'Ερώτημα 2.1α - Διάγραμμα Πόλων-Μηδενικών');
zplane(z, p);
title('Διάγραμμα Πόλων - Μηδενικών (2.1α)');

[b, a] = zp2tf(z, p, K);
disp('[2.1α] Συντελεστές b:'); disp(b);
disp('[2.1α] Συντελεστές a:'); disp(a);

%% 2.1β) Απόκριση Πλάτους και Φάσης
figure('Name', 'Ερώτημα 2.1β - Απόκριση Πλάτους και Φάσης');
freqz(b, a);
title('Απόκριση Συχνότητας Φίλτρου (2.1α)');

%% 2.1γ) Κρουστική και Βηματική Απόκριση
N_samples = 50; 
figure('Name', 'Ερώτημα 2.1γ - Κρουστική και Βηματική Απόκριση');
subplot(2,1,1);
impz(b, a, N_samples);
title('Κρουστική Απόκριση h[n]');
subplot(2,1,2);
stepz(b, a, N_samples);
title('Βηματική Απόκριση s[n]');

%% 2.1δ) Μετακίνηση Πόλων
p_cases = {[0.252 + 0.864i; 0.252 - 0.864i], ...
           [0.280 + 0.960i; 0.280 - 0.960i], ...
           [0.308 + 1.056i; 0.308 - 1.056i]};

figure('Name', 'Ερώτημα 2.1δ - Βηματική Απόκριση και Απόκριση Πλάτους');

for i = 1:3
    p_current = p_cases{i};
    [b_curr, a_curr] = zp2tf(z, p_current, K); 

    
    subplot(3, 2, 2*i - 1);
    stepz(b_curr, a_curr, N_samples);
    title(['Βηματική Απόκριση - Πόλοι ', num2str(i)]);

    
    if i <= 2
        [H, w_f] = freqz(b_curr, a_curr, 512);
        subplot(3, 2, 2*i);
        plot(w_f/pi, abs(H));
        title(['Απόκριση Πλάτους - Πόλοι ', num2str(i)]);
        xlabel('Κανονικοποιημένη Συχνότητα (\times\pi rad/sample)');
        grid on;
    end
end

%% 2.1ε) Διέγερση με Σήμα x[n]
t_sim = 100; 
dt = 1; 
n_vec = 0:t_sim/dt; 
x_n = sin(0.4*pi*n_vec) + sin(0.9*pi*n_vec);
y_n = filter(b, a, x_n); 

figure('Name', 'Ερώτημα 2.1ε - Διέγερση Συστήματος');
plot(n_vec, x_n, 'b', 'DisplayName', 'Σήμα Εισόδου x[n]');
hold on;
plot(n_vec, y_n, 'r', 'DisplayName', 'Σήμα Εξόδου y[n]');
title('Διέγερση Φίλτρου με Σύνθετο Ημιτονοειδές Σήμα');
xlabel('Δείκτης n');
ylabel('Πλάτος');
legend;
grid on;

%% 2.1στ) Εναλλαγή Θέσης Πόλων
p_st = [0.672 + 0.196i; 0.672 - 0.196i];
[b_st, a_st] = zp2tf(z, p_st, K);

figure('Name', 'Ερώτημα 2.1στ - Φίλτρο με Εναλλαγμένους Πόλους');
subplot(2,1,1);
zplane(z, p_st);
title('Διάγραμμα Πόλων - Μηδενικών (2.1στ)');
subplot(2,1,2);
freqz(b_st, a_st);
title('Απόκριση Συχνότητας (2.1στ)');

%% 2.2α & 2.2β) Viola Series
Fs_music = 44100;
try
    [viola_series, Fs_vs] = audioread('viola_series.wav');
catch
     warning('Δεν βρέθηκε το αρχείο viola_series.wav. Παράλειψη 2.2α,β.');
     viola_series = [];
end

if ~isempty(viola_series)
    t_vs = (0:length(viola_series)-1) / Fs_vs;
    figure('Name', 'Ερώτημα 2.2α - Σήμα Viola Series');
    plot(t_vs, viola_series);
    title('Σήμα Viola Series στο Πεδίο του Χρόνου');
    xlabel('Χρόνος (sec)'); ylabel('Πλάτος'); grid on;

    N_fft_vs = Fs_vs; 
    X_vs = fft(viola_series, N_fft_vs); 
    Magnitude_X_vs = abs(X_vs);
    N_half_vs = N_fft_vs/2;
    f_vs = (0:N_half_vs-1) * (Fs_vs/N_fft_vs); 

    figure('Name', 'Ερώτημα 2.2β - Φάσμα Viola Series');
    plot(f_vs, Magnitude_X_vs(1:N_half_vs));
    title('Μέτρο Φάσματος (Γραμμική Κλίμακα) - Viola Series');
    xlabel('Συχνότητα (Hz)'); ylabel('Πλάτος'); grid on;
end

%% 2.2γ) Viola Note και Θεμελιώδης Συχνότητα
try
    [viola_note, Fs_vn] = audioread('viola_note.wav');
catch
    warning('Δεν βρέθηκε το αρχείο viola_note.wav. Παράλειψη 2.2γ,δ.');
    viola_note = [];
end

if ~isempty(viola_note)
    X_vn = fft(viola_note, N_fft_vs); 
    Magnitude_X_vn = abs(X_vn);
    
    figure('Name', 'Ερώτημα 2.2γ - Φάσμα Viola Note');
    plot(f_vs, Magnitude_X_vn(1:N_half_vs));
    title('Μέτρο Φάσματος (Γραμμική Κλίμακα) - Viola Note');
    xlabel('Συχνότητα (Hz)'); ylabel('Πλάτος'); grid on;

    [~, peak_index_vn] = max(Magnitude_X_vn(2:N_half_vs));
    F0_vn = f_vs(peak_index_vn + 1);
    disp(['[2.2γ] Θεμελιώδης Συχνότητα Viola Note: ', num2str(F0_vn), ' Hz']);

    %% 2.2δ) Υλοποίηση Ζωνοπερατών Φίλτρων για Αρμονικές
    F0_example = F0_vn; 
    Fs_vn = Fs_music;
    R = 0.95; 
    
    % 2η Αρμονική
    F2 = 2 * F0_example; 
    omega_c2 = 2*pi * F2 / Fs_vn;
    p_target2 = R * exp(1i * omega_c2); 
    [b2, a2] = zp2tf([0; 0], [p_target2; conj(p_target2)], 1);
    H_2 = freqz(b2, a2, [F2], Fs_vn);
    b2 = b2 .* (1 / abs(H_2(1))); 
    y2 = filter(b2, a2, viola_note);
    
    % 4η Αρμονική
    F4 = 4 * F0_example; 
    omega_c4 = 2*pi * F4 / Fs_vn;
    p_target4 = R * exp(1i * omega_c4); 
    [b4, a4] = zp2tf([0; 0], [p_target4; conj(p_target4)], 1);
    H_4 = freqz(b4, a4, [F4], Fs_vn);
    b4 = b4 .* (1 / abs(H_4(1)));
    y4 = filter(b4, a4, viola_note);
    
    figure('Name', 'Ερώτημα 2.2δ - Απομόνωση 2ης και 4ης Αρμονικής');
    t_seg_vn = (0:length(viola_note)-1)/Fs_vn;
    
    % ΔΙΟΡΘΩΣΗ ΓΙΑ ΤΗ ΓΡΑΜΜΗ 335
    X_y2 = abs(fft(y2, N_fft_vs)); 
    subplot(2, 2, 1);
    plot(f_vs, X_y2(1:N_half_vs));
    title('Φάσμα 2ης Αρμονικής'); xlabel('Hz'); grid on;
    
    subplot(2, 2, 2);
    plot(t_seg_vn(1:1000), y2(1:1000));
    title('Χρόνος 2ης Αρμονικής'); xlabel('sec'); grid on;
    
    % ΔΙΟΡΘΩΣΗ ΓΙΑ ΤΗ ΓΡΑΜΜΗ 341
    X_y4 = abs(fft(y4, N_fft_vs)); 
    subplot(2, 2, 3);
    plot(f_vs, X_y4(1:N_half_vs));
    title('Φάσμα 4ης Αρμονικής'); xlabel('Hz'); grid on;
    
    subplot(2, 2, 4);
    plot(t_seg_vn(1:1000), y4(1:1000));
    title('Χρόνος 4ης Αρμονικής'); xlabel('sec'); grid on;
end
%% ------------------------------------------------------------------------
% ΕΝΟΤΗΤΑ 3: ΔΙΑΧΩΡΙΣΜΟΣ ΜΟΥΣΙΚΩΝ ΝΟΤΩΝ
% ------------------------------------------------------------------------

Fs_mix = 16000;
F0_1 = 370; F0_2 = 440; % F0 for mixture.wav
R = 0.98; % Μέτρο πόλου για στενότερη ζώνη

try
    [mixture, Fs_mix_read] = audioread('mixture.wav');
catch
    warning('Δεν βρέθηκε το αρχείο mixture.wav. Παράλειψη 3α-δ.');
    mixture = [];
end

if ~isempty(mixture)
    %% 3α) Φόρτωση και Φάσμα mixture.wav
    N_fft_mix = Fs_mix; 
    X_mix = fft(mixture, N_fft_mix);
    Magnitude_X_mix = abs(X_mix);
    N_half_mix = N_fft_mix/2;
    f_mix = (0:N_half_mix-1) * (Fs_mix/N_fft_mix);

    figure('Name', 'Ερώτημα 3α - Φάσμα Mixture');
    plot(f_mix, Magnitude_X_mix(1:N_half_mix));
    title('Μέτρο Φάσματος - Mixture.wav');
    xlabel('Συχνότητα (Hz)'); ylabel('Πλάτος'); xlim([0 2000]); grid on;

    %% 3β) Υπολογισμός Συχνοτήτων
    harmonics = 1:4;
    omega_harmonics_1 = 2*pi * (harmonics * F0_1) / Fs_mix;
    omega_harmonics_2 = 2*pi * (harmonics * F0_2) / Fs_mix;
    disp('[3β] Κανονικοποιημένες Συχνότητες (370 Hz):'); disp(omega_harmonics_1);
    disp('[3β] Κανονικοποιημένες Συχνότητες (440 Hz):'); disp(omega_harmonics_2);

    %% 3γ) Σύνθεση Νοτών μέσω Φιλτραρίσματος
    note1_reconstructed = zeros(size(mixture));
    note2_reconstructed = zeros(size(mixture));

    % Νότα 1 (370 Hz)
    for i = 1:length(harmonics)
        omega_c = omega_harmonics_1(i);
        F_h = harmonics(i) * F0_1;
        p_h = R * exp(1i * omega_c);
        [b_h, a_h] = zp2tf([0; 0], [p_h; conj(p_h)], 1);
        H_val = freqz(b_h, a_h, [F_h], Fs_mix);
        b_h = b_h .* (1 / abs(H_val(1)));
        y_h = filter(b_h, a_h, mixture);
        note1_reconstructed = note1_reconstructed + y_h;
    end

    % Νότα 2 (440 Hz)
    for i = 1:length(harmonics)
        omega_c = omega_harmonics_2(i);
        F_h = harmonics(i) * F0_2;
        p_h = R * exp(1i * omega_c);
        [b_h, a_h] = zp2tf([0; 0], [p_h; conj(p_h)], 1);
        H_val = freqz(b_h, a_h, [F_h], Fs_mix);
        b_h = b_h .* (1 / abs(H_val(1)));
        y_h = filter(b_h, a_h, mixture);
        note2_reconstructed = note2_reconstructed + y_h;
    end
    
    note1_reconstructed = note1_reconstructed / max(abs(note1_reconstructed));
    note2_reconstructed = note2_reconstructed / max(abs(note2_reconstructed));


    %% 3δ) Αξιολόγηση Αποτελεσμάτων
    figure('Name', 'Ερώτημα 3δ - Ανακατασκευή Νοτών');
    t_mix = (0:length(mixture)-1)/Fs_mix;
    
    subplot(2, 2, 1);
    plot(t_mix(1:Fs_mix/10), note1_reconstructed(1:Fs_mix/10));
    title('Χρόνος Νότας 1 (370 Hz)'); xlabel('sec'); grid on;
    
    % ΔΙΟΡΘΩΣΗ ΓΙΑ ΤΗ ΓΡΑΜΜΗ 421
    X_note1 = abs(fft(note1_reconstructed, N_fft_mix));
    subplot(2, 2, 3);
    plot(f_mix, X_note1(1:N_half_mix));
    title('Φάσμα Νότας 1'); xlabel('Hz'); xlim([0 2000]); grid on;
    
    subplot(2, 2, 2);
    plot(t_mix(1:Fs_mix/10), note2_reconstructed(1:Fs_mix/10));
    title('Χρόνος Νότας 2 (440 Hz)'); xlabel('sec'); grid on;
    
    % ΔΙΟΡΘΩΣΗ ΓΙΑ ΤΗ ΓΡΑΜΜΗ 428
    X_note2 = abs(fft(note2_reconstructed, N_fft_mix));
    subplot(2, 2, 4);
    plot(f_mix, X_note2(1:N_half_mix));
    title('Φάσμα Νότας 2'); xlabel('Hz'); xlim([0 2000]); grid on;
end
%% 3ε) Επανάληψη για mixture2.wav
F0_3 = 440; F0_4 = 220; % F0 for mixture2.wav

try
    [mixture2, Fs_mix2] = audioread('mixture2.wav');
catch
    warning('Δεν βρέθηκε το αρχείο mixture2.wav.');
    mixture2 = [];
end

if ~isempty(mixture2)
    disp('[3ε] Στην περίπτωση του mixture2.wav, F0,1 = 2*F0,2 (440 Hz = 2*220 Hz).');
    disp('Η απομόνωση των νοτών θα είναι ατελής λόγω της αρμονικής επικάλυψης.');
end