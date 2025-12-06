clear; close all; clc;

%% --- Configuração de figuras proporcionais ---
set(0,'DefaultFigureUnits','pixels');
screen = get(0,'ScreenSize');
defaultFigPos = [125 100 1050 600];

function fh = fixedFig(name, pos)
    scr = get(0,'ScreenSize');
    left   = max(1, min(pos(1), scr(3)-20));
    bottom = max(1, min(pos(2), scr(4)-20));
    width  = min(pos(3), scr(3)-left);
    height = min(pos(4), scr(4)-bottom);
    fh = figure('Name', name, 'NumberTitle','off', ...
                'Units','pixels', 'Position', [left bottom width height]);
    try set(fh,'WindowState','normal'); catch, end
    drawnow;
end

%% ===============================================================
% TRABALHO PRÁTICO 2 - FILTROS FIR E FILTRAGEM ADAPTATIVA
%% ===============================================================

fprintf('========================================\n');
fprintf('TRABALHO PRÁTICO 2\n');
fprintf('Processamento de Sinais\n');
fprintf('========================================\n\n');

%% ===============================================================
% 1. IMPLEMENTAÇÃO BÁSICA
%% ===============================================================

%% ------- 1.1 Carregamento de dados -------
fprintf('=== 1. IMPLEMENTAÇÃO BÁSICA ===\n\n');
fprintf('1.1. Carregando áudio corrompido...\n');

base_dir = fullfile('..','material_fornecido');
audio_file = fullfile(base_dir, 'audio_corrompido.wav');

if ~exist(audio_file, 'file')
    error('Arquivo de audio não encontrado: %s', audio_file);
end

[x, fs] = audioread(audio_file);
if size(x,2) > 1
    x = mean(x,2);  % converte para mono
end
N = length(x);
t = (0:N-1)/fs;

fprintf('   Amostras: %d, fs = %d Hz, duração = %.3f s\n\n', N, fs, N/fs);

% Gráficos do sinal original
fig1 = fixedFig('1. Sinal corrompido - Tempo e Frequência', defaultFigPos);

subplot(3,2,[1 2]);
plot(t, x);
xlabel('Tempo (s)');
ylabel('Amplitude');
title('Forma de onda do áudio corrompido');
grid on;

Nfft = 2^nextpow2(N);
X = fft(x, Nfft);
Xs = fftshift(X);

freqs = (-Nfft/2 : Nfft/2-1) * (fs / Nfft);   % Hz
freqs_khz = freqs / 1000;                     % kHz

amp = 1/N * abs(Xs);
phase = angle(Xs);   % fase embrulhada [–π,π]

% --- Plot amplitude ---
subplot(3,2,3);
plot(freqs_khz, amp, 'b');
xlabel('Frequência (kHz)');
ylabel('|X(e^{j\omega})|');
title('Espectro de amplitude do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);

% --- Plot fase ---
subplot(3,2,4);
plot(freqs_khz, phase, 'b');
xlabel('Frequência (kHz)');
ylabel('Fase (rad)');
title('Espectro de fase do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);

%% ------- EXTRA: Espectrograma do sinal corrompido (Figura 4 do PDF) -------

% Parâmetros do espectrograma (conforme especificação)
window_length_spec = 1024;  % janela de 1024 pontos
overlap_spec = window_length_spec / 2;  % 50% de sobreposição
win_spec = hann(window_length_spec, 'periodic');

% Calcular espectrograma
[S_spec, F_spec, T_spec] = spectrogram(x, win_spec, overlap_spec, window_length_spec, fs);

% Converter para frequência normalizada (0 a 1 em unidades de π rad/amostra) e tempo em amostras
F_norm = F_spec / (fs/2);  % frequência normalizada (0 a 1, onde 1 = π rad/amostra)
T_samples = T_spec * fs;   % tempo em amostras

% Plotar espectrograma
subplot(3,2,[5 6]);
imagesc(T_samples, F_norm, 20*log10(abs(S_spec) + eps));
axis xy;
xlabel('Amostras');
ylabel('Frequência Normalizada (π rad/amostra)');
title('Espectrograma do sinal de áudio corrompido');
colorbar;
colormap('jet');
caxis([-80 0]);
grid on;

%% ------- 1.2 Reprodução do áudio original -------
fprintf('1.2. Reprodução do áudio\n');
choice = input('   Deseja ouvir o áudio corrompido? (s/n): ', 's');
if lower(choice) == 's'
    fprintf('   Reproduzindo áudio corrompido...\n');
    sound(x, fs);
    pause(N/fs + 0.5);
else
    fprintf('   Áudio não reproduzido.\n');
end
fprintf('\n');

%% ------- 1.3 Projeto do Filtro FIR com Janela de Kaiser -------
fprintf('1.3. Projeto do Filtro FIR (Janela de Kaiser)\n');

% Parâmetros do filtro
fp = 5000;  % Hz - frequência de passagem
fr = 6000;  % Hz - frequência de rejeição
Ap = 0.001; % 0.1% = 0.001 em escala linear
Ar = 0.001; % 0.1% = 0.001 em escala linear

% Conversão para dB
delta = min(Ap, Ar);
A = -20*log10(delta);  % Atenuação em dB

fprintf('   fp = %.1f kHz, fr = %.1f kHz\n', fp/1000, fr/1000);
fprintf('   Ap = Ar = %.1f%%\n', Ap*100);
fprintf('   Atenuação A = %.2f dB\n', A);

% Cálculo dos parâmetros da janela de Kaiser
% Oppenheim Seção 7.5.3
if A > 50
    Beta = 0.1102 * (A - 8.7);
elseif A >= 21
    Beta = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21);
else
    Beta = 0.0;
end

% Largura de transição normalizada
delta_f = (fr - fp) / fs;  % normalizada
% Ordem do filtro 
M = ceil((A - 8) / (2.285 * 2 * pi * delta_f));
Num_coef = M + 1;  % M = número de coeficientes

fprintf('   Parâmetro Beta = %.4f\n', Beta);
fprintf('   Ordem do filtro M = %d\n', M);
fprintf('   Número de coeficientes M = %d\n\n', Num_coef);

% Gerar janela de Kaiser (usa M coeficientes)
w_kaiser = kaiser(Num_coef, Beta);

% Filtro ideal passa-baixas
fc = (fp + fr) / 2;  % frequência de corte
wc = 2 * pi * fc / fs;  % frequência angular normalizada

% Resposta ao impulso do filtro ideal (passa-baixas)
n = 0:Num_coef;
alpha = Num_coef/2;  % atraso para tornar causal 
h_ideal = sin(wc * (n - alpha)) ./ (pi * (n - alpha));
h_ideal(n == alpha) = wc / pi;  % corrigir singularidade em n = alpha

% Aplicar janela
h_fir = h_ideal .* w_kaiser';

% Normalizar para ganho unitário em DC
h_fir = h_fir / sum(h_fir);

% Gráfico da janela
fig2 = fixedFig('1.3. Janela de Kaiser', defaultFigPos);
figure(fig2);  % força ativação da janela correta

subplot(2,1,1);
stem(0:N_fir, w_kaiser, 'filled');
xlabel('n');
ylabel('w[n]');
title(sprintf('Janela de Kaiser (M = %d coeficientes, N = %d, β = %.4f)', M, N_fir, Beta));
grid on;

subplot(2,1,2);
stem(0:Num_coef, h_fir, 'filled');
xlabel('n');
ylabel('h[n]');
title('Resposta ao impulso do filtro FIR após janelamento');
grid on;

%% ------- 1.4 Respostas de magnitude e fase do filtro FIR -------
fprintf('1.4. Apresentação das respostas do filtro FIR\n\n');

Nfft_filter = 16384;
[H_fir, w_fir] = freqz(h_fir, 1, Nfft_filter, fs);
w_fir_khz = w_fir / 1000;

fig3 = fixedFig('1.4. Respostas do Filtro FIR', defaultFigPos);
figure(fig3);  % força ativação da janela correta

% Magnitude linear
subplot(2,2,1);
plot(w_fir_khz, abs(H_fir));
xlabel('f (kHz)');
ylabel('|H(f)|');
title('Magnitude (escala linear)');
grid on;
hold on;
xline(fp/1000, 'r--', 'f_P');
xline(fr/1000, 'r--', 'f_R');
hold off;

% Magnitude em dB
subplot(2,2,2);
plot(w_fir_khz, 20*log10(abs(H_fir)));
xlabel('f (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude (escala dB)');
grid on;
ylim([-120 5]);
hold on;
xline(fp/1000, 'r--', 'f_P');
xline(fr/1000, 'r--', 'f_R');
yline(-20*log10(1+Ap), 'k--', sprintf('%.1f dB', -20*log10(1+Ap)));
yline(-A, 'k--', sprintf('%.1f dB', -A));
hold off;

% Fase (unwrapped)
subplot(2,2,3);
plot(w_fir_khz, unwrap(angle(H_fir)));
xlabel('f (kHz)');
ylabel('Fase (rad)');
title('Fase (rad)');
grid on;

% Fase em graus
subplot(2,2,4);
plot(w_fir_khz, unwrap(angle(H_fir))*180/pi);
xlabel('f (kHz)');
ylabel('Fase (graus)');
title('Fase (graus)');
grid on;

%% ------- 1.5 Filtragem usando as 3 abordagens -------
fprintf('1.5. Filtragem do sinal usando 3 abordagens\n');

% Para FIR: num = h_fir, den = 1
num_fir = h_fir;
den_fir = 1;

% (a) Equação de diferenças
fprintf('   a) Filtragem por Equação de Diferenças...\n');
y_eqdif = filtragemPorEqDif(x, num_fir, den_fir);

% (b) Convolução
fprintf('   b) Filtragem por Convolução...\n');
y_conv = conv(x, h_fir);
y_conv = y_conv(1:N);  % truncar ao tamanho original

% (c) FFT
fprintf('   c) Filtragem por FFT...\n');
y_fft = filtragemPorFFT(x, h_fir);
y_fft = y_fft(1:N);  % truncar ao tamanho original

fprintf('   Filtragem concluída.\n\n');

% Visualização
fig4 = fixedFig('1.5. Sinais Filtrados - FIR', defaultFigPos);
figure(fig4);  % força ativação da janela correta

% Eq. Diferenças
subplot(3,3,1);
plot(t, y_eqdif);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Eq. Diferenças (tempo)'); grid on;

Y_eqdif = fftshift(fft(y_eqdif, Nfft))/N;
subplot(3,3,4);
semilogy(freqs_khz, abs(Y_eqdif));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Eq. Diferenças (freq)'); grid on;

subplot(3,3,7);
plot(freqs_khz, angle(Y_eqdif));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Eq. Diferenças (fase)'); grid on;

% Convolução
subplot(3,3,2);
plot(t, y_conv);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Convolução (tempo)'); grid on;

Y_conv = fftshift(fft(y_conv, Nfft))/N;
subplot(3,3,5);
semilogy(freqs_khz, abs(Y_conv));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Convolução (freq)'); grid on;

subplot(3,3,8);
plot(freqs_khz, angle(Y_conv));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Convolução (fase)'); grid on;

% FFT
subplot(3,3,3);
plot(t, y_fft);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('FFT (tempo)'); grid on;

Y_fft = fftshift(fft(y_fft, Nfft))/N;
subplot(3,3,6);
semilogy(freqs_khz, abs(Y_fft));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('FFT (freq)'); grid on;

subplot(3,3,9);
plot(freqs_khz, angle(Y_fft));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('FFT (fase)'); grid on;

%% ------- 1.6 Reprodução do áudio filtrado -------
fprintf('1.6. Reprodução do áudio filtrado\n');
choice = input('   Deseja ouvir o áudio filtrado (FIR)? (s/n): ', 's');
if lower(choice) == 's'
    fprintf('   Reproduzindo áudio filtrado...\n');
    sound(y_eqdif, fs);
    pause(N/fs + 0.5);
else
    fprintf('   Áudio não reproduzido.\n');
end
fprintf('\n');

%% ------- 1.7 Comparação FIR vs IIR -------
fprintf('1.7. Análise Comparativa FIR vs IIR\n');
fprintf('   (Carregar resultados do TP1 para comparação)\n\n');

% Se tiver os coeficientes IIR do TP1:
num_iir_file = fullfile(base_dir, 'coefs_num.mat');
den_iir_file = fullfile(base_dir, 'coefs_den.mat');

if exist(num_iir_file,'file') && exist(den_iir_file,'file')
    fprintf('   Carregando filtro IIR do TP1...\n');
    s_num = load(num_iir_file);
    s_den = load(den_iir_file);
    num_iir = s_num.num;
    den_iir = s_den.den;
    
    % Filtrar com IIR
    y_iir = filtragemPorEqDif(x, num_iir, den_iir);
    
    % Comparação
    fig5 = fixedFig('1.7. Comparação FIR vs IIR', defaultFigPos);
    figure(fig5);  % força ativação da janela correta

    subplot(2,2,1);
    plot(t, y_iir);
    xlabel('Tempo (s)'); ylabel('Amplitude');
    title('Sinal Filtrado - IIR'); grid on;
    
    subplot(2,2,2);
    plot(t, y_eqdif);
    xlabel('Tempo (s)'); ylabel('Amplitude');
    title('Sinal Filtrado - FIR'); grid on;
    
    Y_iir = fftshift(fft(y_iir, Nfft))/N;
    subplot(2,2,3);
    semilogy(freqs_khz, abs(Y_iir));
    hold on;
    semilogy(freqs_khz, abs(Y_eqdif), '--');
    xlabel('f (kHz)'); ylabel('|Y(f)|');
    title('Comparação Espectral');
    legend('IIR', 'FIR');
    grid on;
    
    % Respostas dos filtros
    [H_iir, w_iir] = freqz(num_iir, den_iir, Nfft_filter, fs);
    w_iir_khz = w_iir / 1000;
    
    subplot(2,2,4);
    plot(w_iir_khz, 20*log10(abs(H_iir)));
    hold on;
    plot(w_fir_khz, 20*log10(abs(H_fir)), '--');
    xlabel('f (kHz)'); ylabel('Magnitude (dB)');
    title('Respostas em Frequência');
    legend('IIR', 'FIR');
    grid on;
    ylim([-120 5]);
    
    fprintf('   Comparação visual gerada.\n');
    fprintf('   Observações:\n');
    fprintf('   - FIR: fase linear, sem ripple significativo\n');
    fprintf('   - IIR: ripple ~1dB, fase não-linear, menor ordem\n\n');
else
    fprintf('   Arquivos IIR não encontrados. Comparação não realizada.\n\n');
end

%% ===============================================================
% 2. BÔNUS - FILTRAGEM ADAPTATIVA COM STFT
%% ===============================================================

fprintf('========================================\n');
fprintf('=== 2. BÔNUS - FILTRAGEM ADAPTATIVA ===\n');
fprintf('========================================\n\n');

choice = input('Deseja executar a parte BÔNUS (STFT adaptativa)? (s/n): ', 's');
if lower(choice) ~= 's'
    fprintf('Programa finalizado.\n');
    return;
end

%% ------- 2.1 Cálculo da STFT -------
fprintf('2.1. Cálculo da STFT\n');

% Parâmetros da STFT para reconstrução perfeita
window_length = 1024;  % L
overlap_length = window_length / 2;  % 50% overlap (L/2)
hop_length = window_length - overlap_length;  % R

% CORREÇÃO: FFT maior para melhor resolução
nfft_stft = 2048;  % potência de 2 >= window_length

% Janela de Hann (boa para reconstrução perfeita com 50% overlap)
win = hann(window_length, 'periodic');

fprintf('   Janela: Hann\n');
fprintf('   Comprimento: %d amostras\n', window_length);
fprintf('   Overlap: %d amostras (%.1f%%)\n', overlap_length, 100*overlap_length/window_length);
fprintf('   Hop: %d amostras\n', hop_length);
fprintf('   FFT Length: %d pontos\n\n', nfft_stft);

% Calcular STFT
[S, F, T] = stft(x, fs, 'Window', win, 'OverlapLength', overlap_length, ...
                 'FFTLength', nfft_stft);

fprintf('   Dimensões da STFT: %d x %d (freq x tempo)\n', size(S,1), size(S,2));
fprintf('   Faixa de frequência: 0 a %.1f kHz\n\n', max(F)/1000);

%% ------- 2.2 Espectrograma -------
fprintf('2.2. Apresentação do espectrograma\n\n');

fig6 = fixedFig('2.2. Espectrograma (STFT)', defaultFigPos);
figure(fig6);  % força ativação da janela correta

% 2D
subplot(1,2,1);
imagesc(T, F/1000, 20*log10(abs(S) + eps));
axis xy;
xlabel('Tempo (s)');
ylabel('Frequência (kHz)');
title('Espectrograma 2D');
colorbar;
colormap('jet');
caxis([-80 0]);

% 3D
subplot(1,2,2);
surf(T, F/1000, 20*log10(abs(S) + eps), 'EdgeColor', 'none');
xlabel('Tempo (s)');
ylabel('Frequência (kHz)');
zlabel('Magnitude (dB)');
title('Espectrograma 3D');
view(45, 30);
colorbar;
colormap('jet');
caxis([-80 0]);

%% ------- 2.3 Filtragem Adaptativa -------
fprintf('2.3. Filtragem adaptativa baseada em variância\n');

% Limiar de frequência para análise
f_threshold = fr;  % 6 kHz
idx_threshold = find(F >= f_threshold, 1);

fprintf('   Analisando componentes acima de %.1f kHz\n', f_threshold/1000);
fprintf('   Índice de corte: %d de %d bins de frequência\n', idx_threshold, length(F));

% Analisar variância em cada frame
num_frames = size(S, 2);
variance_high_freq = zeros(1, num_frames);

for k = 1:num_frames
    % Pegar componentes acima do limiar
    high_freq_components = S(idx_threshold:end, k);
    % Calcular variância (energia)
    variance_high_freq(k) = var(abs(high_freq_components));
end

% Detectar frames com ruído (variância alta)
% usar threshold mais robusto
threshold_var = median(variance_high_freq) + 2*std(variance_high_freq);
frames_with_noise = variance_high_freq > threshold_var;

fprintf('   Threshold de variância: %.4e\n', threshold_var);
fprintf('   Frames detectados com ruído: %d de %d (%.1f%%)\n', ...
        sum(frames_with_noise), num_frames, 100*sum(frames_with_noise)/num_frames);

% Calcular resposta do filtro no domínio correto
% A STFT retorna apenas frequências positivas [0, fs/2]
% Precisamos da resposta do filtro nesses pontos
H_fir_stft = freqz(h_fir, 1, length(F), fs);

fprintf('   Dimensões H_fir_stft: %d x 1\n', length(H_fir_stft));
fprintf('   Dimensões S: %d x %d\n\n', size(S,1), size(S,2));

% Aplicar filtro FIR apenas nos frames com ruído
S_filtered = S;  % cópia da STFT

for k = 1:num_frames
    if frames_with_noise(k)
        % multiplicar pela resposta correta
        S_filtered(:, k) = S(:, k) .* H_fir_stft;
    end
end

fprintf('   Filtragem adaptativa concluída.\n\n');

% Visualizar decisão de filtragem
fig7 = fixedFig('2.3. Detecção Adaptativa de Ruído', defaultFigPos);
figure(fig7);  % força ativação da janela correta

subplot(2,1,1);
plot(T, variance_high_freq);
hold on;
yline(threshold_var, 'r--', 'Limiar');
xlabel('Tempo (s)');
ylabel('Variância');
title('Variância das componentes de alta frequência');
grid on;
legend('Variância', 'Threshold');

subplot(2,1,2);
imagesc(T, F/1000, 20*log10(abs(S_filtered) + eps));
axis xy;
xlabel('Tempo (s)');
ylabel('Frequência (kHz)');
title('Espectrograma após filtragem adaptativa');
colorbar;
colormap('jet');
caxis([-80 0]);
hold on;
% Marcar frames filtrados
for k = 1:num_frames
    if frames_with_noise(k)
        plot([T(k) T(k)], [0 max(F)/1000], 'w--', 'LineWidth', 0.5);
    end
end

%% ------- 2.4 Reconstrução (iSTFT) -------
fprintf('2.4. Reconstrução do sinal (iSTFT)\n\n');

y_adaptive = istft(S_filtered, fs, 'Window', win, 'OverlapLength', overlap_length, ...
                   'FFTLength', nfft_stft);

% Ajustar tamanho
if length(y_adaptive) > N
    y_adaptive = y_adaptive(1:N);
elseif length(y_adaptive) < N
    y_adaptive = [y_adaptive; zeros(N - length(y_adaptive), 1)];
end

fprintf('   Sinal reconstruído: %d amostras\n\n', length(y_adaptive));

% Visualização
fig8 = fixedFig('2.4. Sinal Reconstruído - Filtragem Adaptativa', defaultFigPos);
figure(fig8);  % força ativação da janela correta

subplot(3,1,1);
plot(t, x);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Sinal Original (corrompido)');
grid on;

subplot(3,1,2);
plot(t, y_eqdif);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Filtragem Linear (FIR)');
grid on;

subplot(3,1,3);
disp(max(abs(imag(y_adaptive))))
plot(t, y_adaptive);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Filtragem Adaptativa (STFT)');
grid on;

%% ------- 2.5 Reprodução do áudio adaptativo -------
fprintf('2.5. Reprodução do áudio com filtragem adaptativa\n');
choice = input('   Deseja ouvir o áudio filtrado adaptativamente? (s/n): ', 's');
if lower(choice) == 's'
    fprintf('   Reproduzindo áudio adaptativo...\n');
    sound(y_adaptive, fs);
    pause(N/fs + 0.5);
else
    fprintf('   Áudio não reproduzido.\n');
end
fprintf('\n');

%% ------- 2.6 Análise Comparativa -------
fprintf('2.6. Análise Comparativa Final\n\n');

fig9 = fixedFig('2.6. Comparação Final: Linear vs Adaptativa', defaultFigPos);
figure(fig9);  % força ativação da janela correta

% Espectros
Y_linear = fftshift(fft(y_eqdif, Nfft))/N;
Y_adaptive = fftshift(fft(y_adaptive, Nfft))/N;

subplot(2,2,1);
semilogy(freqs_khz, abs(Y_linear));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Espectro - Filtragem Linear (FIR)');
grid on;
xlim([0 10]);

subplot(2,2,2);
semilogy(freqs_khz, abs(Y_adaptive));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Espectro - Filtragem Adaptativa');
grid on;
xlim([0 10]);

subplot(2,2,3);
plot(freqs_khz, abs(Y_linear));
hold on;
plot(freqs_khz, abs(Y_adaptive), '--');
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Comparação Espectral');
legend('Linear', 'Adaptativa');
grid on;
xlim([0 10]);

subplot(2,2,4);
plot(t, abs(y_eqdif - y_adaptive));
xlabel('Tempo (s)'); ylabel('Diferença Absoluta');
title('Diferença entre Filtragem Linear e Adaptativa');
grid on;

% Calcular métricas
snr_linear = 10*log10(sum(y_eqdif.^2) / sum((x-y_eqdif).^2));
snr_adaptive = 10*log10(sum(y_adaptive.^2) / sum((x-y_adaptive).^2));

fprintf('Métricas de Desempenho:\n');
fprintf('   SNR Filtragem Linear:    %.2f dB\n', snr_linear);
fprintf('   SNR Filtragem Adaptativa: %.2f dB\n', snr_adaptive);
fprintf('   Diferença:                %.2f dB\n\n', snr_adaptive - snr_linear);

fprintf('Análise:\n');
fprintf('   - Filtragem linear: remove ruído globalmente\n');
fprintf('   - Filtragem adaptativa: preserva alta frequência fora do ruído\n');
fprintf('   - Comparar espectros e audição para avaliar qualidade\n\n');

%% ===============================================================
% FUNÇÕES AUXILIARES - TRABALHO PRÁTICO 1
%% ===============================================================

function y = filtragemPorEqDif(x, num, den)
    % Filtragem por equação de diferenças
    % Normalização para garantir que den(1) = 1
    if den(1) ~= 1
        num = num / den(1);
        den = den / den(1);
    end
    
    Nx = length(x);
    M = length(num);
    N = length(den);
    y = zeros(size(x));
    
    % Laço amostra a amostra
    for n = 1:Nx
        acc = 0;
        
        % Soma dos termos num(k)*x[n-k+1]
        for k = 1:M
            idx = n - k + 1;
            if idx > 0
                acc = acc + num(k) * x(idx);
            end
        end
        
        % Feedback: subtração dos termos den(k)*y[n-k+1], k>=2
        for k = 2:N
            idx = n - k + 1;
            if idx > 0
                acc = acc - den(k) * y(idx);
            end
        end
        
        y(n) = acc;
    end
end

function y = filtragemPorFFT(x, h)
    % Filtragem pela multiplicação da FFT
    % Comprimentos
    Nx = length(x);
    Nh = length(h);
    
    % Tamanho da FFT (potência de 2 >= Nx+Nh-1)
    % Limitar tamanho máximo para evitar overflow
    Nfft_ideal = Nx + Nh - 1;
    Nfft = 2^(nextpow2(Nfft_ideal));
    
    % Verificar se não excede limite
    max_size = 2^27;  % ~134 milhões de elementos
    if Nfft > max_size
        fprintf('   AVISO: FFT muito grande (%d), usando tamanho limitado\n', Nfft);
        Nfft = max_size;
    end
    
    % Garantir que x e h são vetores coluna
    x = x(:);
    h = h(:);
    
    % FFT do sinal e da resposta
    X = fft(x, Nfft);
    H = fft(h, Nfft);
    
    % Multiplicação no domínio da frequência
    Y = X .* H;
    
    % IFFT e truncagem para tamanho correto
    y = real(ifft(Y));
    
    % Retornar apenas a parte válida da convolução
    y = y(1:min(Nfft_ideal, length(y)));
    
    % Garantir formato de coluna
    y = y(:);
end