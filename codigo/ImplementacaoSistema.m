clear; close all; clc;

%% --- Forçar figuras proporcionais em qualquer máquina ---
set(0,'DefaultFigureUnits','pixels');   % força uso de pixels
screen = get(0,'ScreenSize');          % [left bottom width height] da tela

% Posição/tamanho 
defaultFigPos = [125 100 1050 600];   % [left bottom width height]

% Função auxiliar para criar figuras com posição e tamanho fixos
function fh = fixedFig(name, pos)
    scr = get(0,'ScreenSize');
    % garante que a janela caiba na tela 
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
% 1. Carregamento de dados e apresentação de características básicas 
%% ===============================================================
%% ------- Parâmetros -------
base_dir = fullfile('..','material_fornecido');

audio_file = fullfile(base_dir, 'audio_corrompido.wav');
num_file   = fullfile(base_dir, 'coefs_num.mat');
den_file   = fullfile(base_dir, 'coefs_den.mat');
nfft_spec  = 16384;    % pontos para FFT. Poderia ser outra potência de 2, mas 16384 é um compromisso:
%suficientemente grande para ter boa resolução em Hz, mas sem deixar a FFT lenta

%% ------- 1.1 Carregamento do áudio e reprodução -------
if ~exist(audio_file, 'file')
    error('Arquivo de audio não encontrado: %s', audio_file);
end

[x, fs] = audioread(audio_file);      % lê áudio
if size(x,2) > 1
    x = mean(x,2);                    % converte para mono (média dos canais)
end
N = length(x);
t = (0:N-1)/fs;

%Extra apenas
%fprintf('Amostras: %d, fs = %d Hz, duração = %.3f s\n', N, fs, N/fs);

%% Perguntar ao usuário se deseja ouvir
fprintf('TRABALHO PRÁTICO 1\n=== Seção 1: Carregamento de dados e apresentação de características básicas  ===\n');

choice = input('Deseja ouvir o áudio original (corrompido)? (s/n): ', 's');

if lower(choice) == 's'
    try
        fprintf('Reproduzindo áudio original (corrompido)... (aguarde)\n');
        sound(x, fs);
        pause(N/fs); % toca até o fim
    catch ME
        warning(ME.identifier, 'Não foi possível reproduzir áudio: %s', ME.message);
    end
else
    fprintf('Ok, o áudio não será reproduzido.\n');
end

%% ------- 1.2 Gráfico da forma de onda (tempo em s) -------
fig1 = fixedFig('1. Carregamento de dados e apresentação de características básicas', defaultFigPos);
figure(fig1);  % força ativação da janela correta

subplot(5,2,[1 2]); % pega colunas 1 e 2 da primeira linha
plot(t, x);
xlabel('Tempo (s)');
ylabel('Amplitude normalizada');
title('Forma de onda do áudio corrompido');
grid on;

%% ------- 1.3 Espectros de amplitude |X(e^{j\omega})|×f (kHz) e fase θ(ω)×f (kHz) -------
Nfft = max(nfft_spec, 2^nextpow2(N));
X = fft(x, Nfft);
Xs = fftshift(X);

freqs = (-Nfft/2 : Nfft/2-1) * (fs / Nfft);   % Hz
freqs_khz = freqs / 1000;                     % kHz

amp = 1/N * abs(Xs);
phase = angle(Xs);   % fase embrulhada [–π,π]

% --- Plot amplitude ---
subplot(5,2,3);
plot(freqs_khz, amp, 'b');
xlabel('Frequência (kHz)');
ylabel('|X(e^{j\omega})|');
title('Espectro de amplitude do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);

% --- Plot fase ---
subplot(5,2,4);
plot(freqs_khz, phase, 'b');
xlabel('Frequência (kHz)');
ylabel('Fase (rad)');
title('Espectro de fase do áudio corrompido');
grid on;
xlim([min(freqs_khz) max(freqs_khz)]);


%% ------- 1.4 Carregamento dos coeficientes do filtro -------
if ~exist(num_file,'file') || ~exist(den_file,'file')
    error('Arquivos de coeficientes não encontrados. Esperados: %s e %s', num_file, den_file);
end
s_num = load(num_file);    % coeficientes numerador
s_den = load(den_file);    % coeficientes denominador 

num = s_num.num; % conteúdo do campo chamado num dentro da struct s_num.
den = s_den.den;

%% EXTRA: Impressão da função de transferência H(z)
%fprintf('\n=== Função de transferência H(z) ===\n');
%Hz = tf(num, den, -1)   % -1 força o domínio z^-1

% Avaliar no ponto z = 1
%z = 1;
%H1 = polyval(num, z) / polyval(den, z);

%fprintf('H(z=1) = %.6e\n', H1);
% ver tamanho/estrutura
%whos num den

%% ------- 1.5 Respostas de magnitude e fase |H(e^{jω})| e θ(ω) -------

Nfft = 16384;

% ----------- a) Escala linear simétrica (-fs/2..+fs/2) ----------
[H_whole, w_whole] = freqz(num, den, Nfft, 'whole', fs); % 'whole': 0 a fs
Hshift = fftshift(H_whole);  % -fs/2 a fs/2
wshift = (-fs/2 : fs/Nfft : fs/2 - fs/Nfft);
wshift_khz = wshift/1000;

% Magnitude linear
subplot(5,2,5);
plot(wshift_khz, abs(Hshift));
xlabel('f (kHz)');
ylabel('|H(e^{j\omega})|');
title('Magnitude (linear)');
grid on;

% Fase linear
subplot(5,2,7);
plot(wshift_khz, unwrap(angle(Hshift))); % fase desembrulhada 
xlabel('f (kHz)');
ylabel('\theta(\omega) (rad)');          
title('Fase');
grid on;

% ----------- b) Escala log (dB) 0..fs/2 ----------
[H, w] = freqz(num, den, Nfft, fs);  % 0 a fs/2.
w_khz = w/1000;

% Magnitude em dB
subplot(5,2,6);
plot(w_khz, 20*log10(abs(H)));
xlabel('Frequência (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude (dB)');
grid on;
ylim([-420 20]); % para destacar

% Fase em graus (unwrap)
subplot(5,2,8);
plot(w_khz, unwrap(angle(H))*180/pi);
xlabel('Frequência (kHz)');
ylabel('Fase (graus)');
title('Fase');
grid on;

%% ------- 1.6 Resposta ao impulso do filtro -------
n_samples = 1000;   % quantidade de amostras da resposta ao impulso
[h, n] = impz(num, den, n_samples);

% Plotar a resposta ao impulso
subplot(5,2,[9 10]);
stem(n, h, 'filled');
xlabel('n');
ylabel('h[n]');
title('Resposta ao impulso do filtro');
grid on;


%% ===============================================================
% 2. Implementação de funções para filtragem
%% ===============================================================
%% ------- 2.1 Filtragem pela equação de diferenças -------
function y = filtragemPorEqDif(x, num, den)
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
        acc = 0;              % acumulador para calcular y[n]

        % --- soma dos termos num(k)*x[n-k+1] ---
        for k = 1:M
            idx = n - k + 1;         % índice correspondente em x (MATLAB 1-based)
            if idx > 0               % somente se índice válido (condição inicial = 0)
                acc = acc + num(k) * x(idx);
            end
        end

        % --- Feedback: subtração dos termos den(k)*y[n-k+1], k>=2 ---
        % (começa em k=2 porque den(1), depois da normalização certamente unitário, corresponde ao coeficiente de y[n])
        for k = 2:N
            idx = n - k + 1;         % índice correspondente em y
            if idx > 0               % somente se índice válido (condição inicial = 0)
                acc = acc - den(k) * y(idx);
            end
        end

        % --- Armazena resultado ---
        y(n) = acc;
    end
end

%% ------- 2.2 Filtragem pela convolução com resposta ao impulso -------
%% (a) Truncagem da resposta ao impulso
function h_trunc = truncarResposta(h)
    % Critério: manter todas as amostras até o último índice cujo valor
    % seja maior ou igual a 1% do valor de pico da resposta.
    limiar = 0.01 * max(abs(h));
    idx = find(abs(h) >= limiar, 1, 'last');
    h_trunc = h(1:idx);
end

h_trunc = truncarResposta(h);

%% (b) Apresentação da resposta truncada
fig2 = fixedFig('2. Resposta ao impulso truncada', defaultFigPos);
figure(fig2);  % força ativação da janela correta

Nh = length(h_trunc); % Número de amostras após truncagem: Nh = 277
fprintf('\n=== Seção 2: Implementação de funções para filtragem de 3 formas distintas ===\n');
fprintf('Subseção 2.2(b)\n');
fprintf('Resposta truncada: Nh = %d amostras\n', Nh);

stem(0:Nh-1, h_trunc, 'filled');
xlabel('n');
ylabel('h_{trunc}[n]');
title(sprintf('Resposta ao impulso truncada (Nh = %d)', Nh));
grid on;

%% (c) Filtragem por convolução circular
function y = filtragemPorConv(x, h)
    % A convolução circular é feita pelo comando cconv
    % O tamanho é definido como Nx + Nh - 1
    Nx = length(x);
    Nh = length(h);
    N  = Nx + Nh - 1;  % mesmo tamanho pedido no enunciado

    y = cconv(x, h, N);
    y = y(:);  % garante vetor coluna
    
    % Obs.: justo usar cconv para efeitos de comparação, uma vez que 
    % a filtragem por FFT também se utiliza de funções nativas (fft e ifft)

    % De qualquer forma, segue o Método "bruto" sem as otimizações da função cconv:
    % x = x(:).';
    % h = h(:).';
    % 
    % Nx = length(x);
    % Nh = length(h);
    % Ny = Nx + Nh - 1;
    % 
    % % Zero-padding até Ny
    % x_pad = [x zeros(1, Ny - Nx)];
    % h_pad = [h zeros(1, Ny - Nh)];
    % 
    % % Convolução circular
    % y = zeros(1, Ny);
    % for n = 0:Ny-1
    %     y(n+1) = sum(x_pad .* circshift(fliplr(h_pad), [0 n]));
    % end
    % 
    % % Garante saída como vetor coluna
    % y = y(:);
end

%% ------- 2.3 Filtragem pela multiplicação da FFT -------
function y = filtragemPorFFT(x, h)
    % Comprimentos
    Nx = length(x);
    Nh = length(h);

    % Tamanho da FFT (potência de 2 >= Nx+Nh-1)
    Nfft = 2^(nextpow2(Nx+Nh-1));  

    % FFT do sinal e da resposta truncada
    X = fft(x, Nfft);
    H = fft(h, Nfft);

    % Multiplicação no domínio da frequência
    Y = X .* H;

    % IFFT e truncagem para tamanho correto
    y = (real(ifft(Y)));
    y = y(1:Nx+Nh-1);
end


%% ===============================================================
% 3. Filtragem do sinal
%% ===============================================================
fprintf('\n=== Seção 3: Filtragem do sinal  ===\n');

% --- 3.1 Filtragem pelas três formas ---

% Eq. Diferenças
y_eqdif = filtragemPorEqDif(x, num, den);

% Convolução com resposta truncada
y_conv = filtragemPorConv(x, h_trunc);

% FFT
y_fft = filtragemPorFFT(x, h_trunc);

%% --- 3.2 Apresentação dos sinais filtrados no tempo e na frequência ---
fig3 = fixedFig('3. Filtragem do sinal', defaultFigPos);
figure(fig3);  % força ativação da janela correta

% --- Eq. Diferenças ---
subplot(3,3,1);
plot(t, y_eqdif);
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Eq. Diferenças (tempo)'); grid on;

freqs_eq = linspace(-fs/2, fs/2, length(y_eqdif));
freqs_khz_eq = freqs_eq / 1000;
subplot(3,3,4);
Yeq = fftshift(fft(y_eqdif))/length(y_eqdif); 
semilogy(freqs_khz_eq, abs(Yeq));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Eq. Diferenças (freq)'); grid on;

subplot(3,3,7);
plot(freqs_khz, angle(Yeq));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Eq. Diferenças (fase)'); grid on;

% --- Convolução ---
subplot(3,3,2);
plot(t, y_conv(1:N)); 
xlabel('Tempo (s)'); ylabel('Amplitude');
title('Convolução (tempo)'); grid on;

freqs_c = linspace(-fs/2, fs/2, length(y_conv));
freqs_khz_c = freqs_c / 1000;
subplot(3,3,5);
Yc = fftshift(fft(y_conv))/length(y_conv); 
semilogy(freqs_khz_c, abs(Yc));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('Convolução (freq)'); grid on;

subplot(3,3,8);
plot(freqs_khz_c, angle(Yc));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('Convolução (fase)'); grid on;

% --- FFT ---
subplot(3,3,3);
plot(t, y_fft(1:N));
xlabel('Tempo (s)'); ylabel('Amplitude');
title('FFT (tempo)'); grid on;

freqs_fft = linspace(-fs/2, fs/2, length(y_fft));
freqs_khz_fft = freqs_fft / 1000;
subplot(3,3,6);
Yf = fftshift(fft(y_fft))/length(y_fft); 
semilogy(freqs_khz_fft, abs(Yf));
xlabel('f (kHz)'); ylabel('|Y(f)|');
title('FFT (freq)'); grid on;

subplot(3,3,9);
plot(freqs_khz_fft, angle(Yf));
xlabel('f (kHz)'); ylabel('Fase (rad)');
title('FFT (fase)'); grid on;

%% --- 3.3 Reprodução dos sinais filtrados ---

tocou = false; % flag para saber se algum áudio foi reproduzido

choice = input('a) Deseja ouvir o áudio após filtragem por Equação de Diferenças? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (Equação de Diferenças)...\n');
    sound(y_eqdif, fs); pause(N/fs);
    tocou = true;
end

choice = input('b) Deseja ouvir o áudio após filtragem por Convolução? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (Convolução)...\n');
    sound(y_conv, fs); pause(N/fs);
    tocou = true;
end

choice = input('c) Deseja ouvir o áudio após filtragem por FFT? (s/n): ','s');
if lower(choice) == 's'
    fprintf('Reproduzindo áudio filtrado (FFT)...\n');
    sound(y_fft, fs); pause(N/fs);
    tocou = true;
end

if ~tocou
    fprintf('Ok, nenhum áudio será reproduzido.\n');
end


%% ===============================================================
% 4. BÔNUS: Implementação overlap-add 
%         + Benchmarks (EqDif, Conv, FFT, OL-Conv, OL-FFT)
%% ===============================================================
fprintf('\n=== Seção 4: Bônus - Overlap-Add (Convolução e FFT)  ===\n');  % imprime cabeçalho indicando início da seção 4
% Funções auxiliares: Overlap-Add Conv e FFT
function y_out = filtragemPorConvOL(x, h_trunc, L, Ny)
    % ----------------- Overlap-Add usando convolução circular -----------------
    Nx = length(x);
    y_out = zeros(Ny,1);  % inicializa vetor de saída com zeros

    for k = 0:ceil(Nx/L)-1
        idx0 = k*L + 1;                  % índice inicial do bloco
        idx1 = min((k+1)*L, Nx);         % índice final (não ultrapassar Nx)
        xblk = x(idx0:idx1);             % bloco atual de entrada

        % ---------- Convolução circular (mesmo código de filtragemPorConv) ----------
        Nx_blk = length(xblk);
        Nh = length(h_trunc);
        N  = Nx_blk + Nh - 1;            % tamanho da convolução circular
        yblk = cconv(xblk, h_trunc, N);  % filtragem circular
        yblk = yblk(:);                  % vetor coluna
        % ---------------------------------------------------------------------------

        % ---------- Etapa de Overlap-Add ----------
        out_start = idx0;                                % posição inicial na saída
        out_end = min(Ny, idx0 + length(yblk) - 1);       % posição final (sem ultrapassar Ny)
        len_write = out_end - out_start + 1;              % nº de amostras a escrever
        y_out(out_start:out_end) = y_out(out_start:out_end) + yblk(1:len_write);
    end
end


function y_out = filtragemPorFFTOl(x, h_trunc, L, Ny, Nfft_blk, Hfft_fixed)
    % ----------------- Overlap-add FFT -----------------
    Nx = length(x);
    Nh = length(h_trunc);
    y_out = zeros(Ny,1);                                % inicializa saída acumulada
    for k = 0:ceil(Nx/L)-1                              % para cada bloco, L é o tamanho do bloco
        idx0 = k*L + 1;                                 % índice inicial (1-based)
        idx1 = min((k+1)*L, Nx);                        % índice final do bloco
        xblk = x(idx0:idx1);                            % extrai bloco de x
        Nblk = length(xblk);                            % comprimento do bloco (pode ser < L no último)
        Xblk = fft(xblk, Nfft_blk);                     % FFT do bloco (com zero-padding até Nfft_blk)
        Yblk = ifft(Xblk .* Hfft_fixed);                % multiplicação no domínio freq e IFFT
        yblk = real(Yblk(1:Nblk+Nh-1));                 % extrai a parte linear relevante e força real
        out_start = idx0;                               % posição inicial de escrita
        out_end = min(Ny, idx0 + length(yblk) - 1);     % posição final de escrita (limitada por Ny)
        len_write = out_end - out_start + 1;            % quantas amostras escrever
        y_out(out_start:out_end) = y_out(out_start:out_end) + yblk(1:len_write); % soma (overlap-add)
    end
end

% ----------------- Preparações -----------------
% garantir h_trunc (se não existir, gera)
if ~exist('h_trunc','var') || isempty(h_trunc)            % verifica se a variável h_trunc existe e não está vazia
    fprintf('h_trunc não encontrado — gerando...\n'); % informa que vai gerar h_trunc 
    h_trunc = truncarResposta(h); % chama a função que gera a resposta truncada e guarda em h_trunc
end

Nx = length(x);                   % comprimento do sinal de entrada x
Nh = length(h_trunc);             % comprimento da resposta truncada h_trunc
if Nh <= 0                        % checagem: Nh deve ser positivo
    error('h_trunc vazio ou Nh <= 0. Verifique truncagem em 2.2.'); % lança erro se inválido
end
Ny = Nx + Nh - 1;                 % comprimento resultante da convolução (saída completa)

% Parâmetros de benchmark
reps = 6;    % número de repetições para média de tempo

% --- Warm-up (uma execução de cada para JIT/cache) ---
y_tmp1 = filtragemPorEqDif(x, num, den);              % execução única de EqDif para "aquecer" JIT/caches
y_tmp2 = filtragemPorConv(x, h_trunc);          % execução única de Conv para warm-up
y_tmp3 = filtragemPorFFT(x, h_trunc);                % execução única de FFT para warm-up

% ----------------- Benchmarks principais -----------------
% EqDif
tic;                                                 % inicia cronômetro
for k = 1:reps                                       % repete o processamento reps vezes
    y_eq_bench = filtragemPorEqDif(x, num, den);     % executa filtragem por equação de diferenças
end
t_eq = toc / reps;                                   % tempo médio por execução (s)

% Conv (usa h_trunc)
tic;                                                 % inicia cronômetro
for k = 1:reps                                       % repete reps vezes
    y_conv_bench = filtragemPorConv(x, h_trunc);     % executa filtragem por convolução (usa h_trunc)
end
t_conv = toc / reps;                                 % tempo médio por execução (s)

% FFT (multiplicação FFT)
tic;                                                 % inicia cronômetro
for k = 1:reps                                       % repete reps vezes
    y_fft_bench = filtragemPorFFT(x, h_trunc);       % executa filtragem via FFT (multiplicação no domínio freq)
end
t_fft = toc / reps;                                  % tempo médio por execução (s)

% ----------------- Overlap-Add (seguindo enunciado L = Nh) -----------------
L = Nh;                      % tamanho do bloco de entrada L (aqui igual a Nh)
numBlocks = ceil(Nx / L);     % número de blocos para cobrir todo x     
Nfft_blk = 2^(nextpow2(L + Nh - 1));   % tamanho da FFT por bloco (potência de 2)
Hfft_fixed = fft(h_trunc, Nfft_blk);  % FFT de h_trunc, calculada uma vez 
% Calcular Hfft_fixed aqui e não em filtragemPorFFTOl está ajudando
% o tempo de execução de filtragemPorFFTOl, assim como passar h_trunc como
% argumento já ajuda as outras funções também

% Overlap-add Convolução benchmark
tic;                                                      % inicia cronômetro
for r = 1:reps                                            % repete reps vezes para média
    y_ol_conv_bench = filtragemPorConvOL(x, h_trunc, L, Ny); % chamada da função auxiliar
end
t_ol_conv = toc / reps;                                   % tempo médio OL-Conv (s)

%  Overlap-add FFT benchmark
tic;                                                      % inicia cronômetro
for r = 1:reps                                            % repete reps vezes
    y_ol_fft_bench = filtragemPorFFTOl(x, h_trunc, L, Ny, Nfft_blk, Hfft_fixed); % chamada da função auxiliar
end
t_ol_fft = toc / reps;                                    % tempo médio OL-FFT (s)

% ----------------- Resultados dos benchmarks -----------------
fprintf('Tempos médios (s):\n  Equação de Diferenças = %.6f\n  Convolução  = %.6f\n  FFT   = %.6f\n  Convolução Overlap-add = %.6f\n  FFT  Overlap-add = %.6f\n', ...
    t_eq, t_conv, t_fft, t_ol_conv, t_ol_fft);            % imprime os tempos médios formatados

% ----------------- Gerar saídas "reais" (uma execução final) -----------------
% executar métodos uma vez para comparar sinais/resultados
y_eqdif = filtragemPorEqDif(x, num, den);                 % versão final por equação de diferenças
y_conv = filtragemPorConv(x, h_trunc);               % convolução (também atualiza h_trunc caso necessário)
y_fft = filtragemPorFFT(x, h_trunc);                      % filtragem por FFT (versão final)

% Executar overlap-add (uma vez) para obter y_ol_conv e y_ol_fft finais
y_ol_conv = filtragemPorConvOL(x, h_trunc, L, Ny);        % execução final OL-Conv
y_ol_fft  = filtragemPorFFTOl(x, h_trunc, L, Ny, Nfft_blk, Hfft_fixed); % execução final OL-FFT

% ----------------- Comparações numéricas -----------------
% garantir mesmas dimensões
y_eqdif_tr = [y_eqdif(:); zeros(max(0, Ny - length(y_eqdif)), 1)];  % pad com zeros se necessário e garante coluna
y_eqdif_tr = y_eqdif_tr(1:Ny);                                      % trunca/corta para comprimento Ny

y_conv_tr = [y_conv(:); zeros(max(0, Ny - length(y_conv)), 1)];     % mesma pad/truncagem para y_conv
y_conv_tr = y_conv_tr(1:Ny);                                        % garante Ny

y_fft_tr  = [y_fft(:); zeros(max(0, Ny - length(y_fft)), 1)];       % pad/truncagem para y_fft
y_fft_tr  = y_fft_tr(1:Ny);                                         % garante Ny

y_ol_conv = y_ol_conv(:);                                           % garante vetor coluna
y_ol_fft  = y_ol_fft(:);                                            % garante vetor coluna

% Comparações (erros)
maxabs_olfft_eqdif = max(abs(y_ol_fft - y_eqdif_tr));                
rmsrel_olfft_eqdif = norm(y_ol_fft - y_eqdif_tr) / max(eps, norm(y_eqdif_tr)); 

maxabs_olconv_eqdif = max(abs(y_ol_conv - y_eqdif_tr));              
rmsrel_olconv_eqdif = norm(y_ol_conv - y_eqdif_tr) / max(eps, norm(y_eqdif_tr)); 

maxabs_olfft_fft = max(abs(y_ol_fft - y_fft_tr));                
rmsrel_olfft_fft = norm(y_ol_fft - y_fft_tr) / max(eps, norm(y_fft_tr)); 

maxabs_olconv_conv = max(abs(y_ol_conv - y_conv_tr));              
rmsrel_olconv_conv = norm(y_ol_conv - y_conv_tr) / max(eps, norm(y_conv_tr)); 

fprintf('\nComparação numérica - Overlap-add vs Equação de Diferenças:\n');           
fprintf('Overlap-FFT:  erro máximo absoluto = %.4e, erro RMS relativo = %.4e\n', maxabs_olfft_eqdif, rmsrel_olfft_eqdif); 
fprintf('Overlap-Conv: erro máximo absoluto = %.4e, erro RMS relativo = %.4e\n', maxabs_olconv_eqdif, rmsrel_olconv_eqdif); 
fprintf('\nComparação numérica - Overlap-add vs métodos diretos:\n');           
fprintf('Overlap-FFT vs FFT direta:  erro máximo absoluto = %.4e, erro RMS relativo = %.4e\n', maxabs_olfft_fft, rmsrel_olfft_fft); 
fprintf('Overlap-Conv vs Conv direta: erro máximo absoluto = %.4e, erro RMS relativo = %.4e\n\n', maxabs_olconv_conv, rmsrel_olconv_conv); 

% ----------------- Plots finais -----------------
Y_olfft = fftshift(fft(y_ol_fft)) / length(y_ol_fft);     % espectro OL-FFT
Y_olconv = fftshift(fft(y_ol_conv)) / length(y_ol_conv);  % espectro OL-Conv

% Eixos de frequência baseados no tamanho real de cada sinal
freqs_olfft = linspace(-fs/2, fs/2, length(y_ol_fft));
freqs_khz_olfft = freqs_olfft / 1000;

freqs_olconv = linspace(-fs/2, fs/2, length(y_ol_conv));
freqs_khz_olconv = freqs_olconv / 1000;

fig4 = fixedFig('4. BÔNUS: Implementação overlap-add', defaultFigPos);
figure(fig4);  % força ativação da janela correta

tiledlayout(3,3, 'Padding','compact', 'TileSpacing','compact'); % cria figura para comparação

nexttile;                                                    % posição 1ª linha, 1ª coluna
plot(t, y_ol_conv(1:N));                                           % plota tempo da saída OL-Conv (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Amplitude'); title('Overlap-Add (conv) - tempo'); grid on;

nexttile;                                                  % posição 1ª linha, 2ª coluna 
plot(t, y_ol_fft(1:N));                                            % plota tempo da saída OL-FFT (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Amplitude'); title('Overlap-Add (FFT) - tempo'); grid on;

nexttile;                                                  % posição 1ª linha, 3ª coluna 
plot(t, (y_ol_conv(1:N) - y_eqdif_tr(1:N)));                       % plota erro no tempo entre OL-Conv e EqDif (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Erro'); title('Erro Overlap-Conv - EqDif'); grid on;

nexttile;                                                  % posição 2ª linha, 1ª coluna
semilogy(freqs_khz_olconv, abs(Y_olconv));                      % plota magnitude do espectro OL-Conv
xlabel('f (kHz)'); ylabel('|Y|'); title('Overlap-Add (conv) - freq'); grid on;

nexttile;                                                  % posição 2ª linha, 2ª coluna
semilogy(freqs_khz_olfft, abs(Y_olfft));                       % plota magnitude do espectro OL-FFT
xlabel('f (kHz)'); ylabel('|Y|'); title('Overlap-Add (FFT) - freq'); grid on;

nexttile;                                                  % posição 2ª linha, 3ª coluna 
plot(t, (y_ol_fft(1:N) - y_eqdif_tr(1:N)));                        % plota erro no tempo entre OL-FFT e EqDif (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Erro'); title('Erro Overlap-FFT - EqDif '); grid on;

nexttile;                                                   % posição 3ª linha, 1ª coluna 
plot(t, (y_ol_conv(1:N) - y_conv_tr(1:N)));                        % plota erro entre OL-Conv e convolução direta (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Erro'); title('Erro Overlap-Conv - Conv direta'); grid on;

nexttile;                                                   % posição 3ª linha, 2ª coluna 
plot(t, (y_ol_fft(1:N) - y_fft_tr(1:N)));                          % plota erro entre OL-FFT e FFT direta (apenas N amostras)
xlabel('Tempo (s)'); ylabel('Erro'); title('Erro Overlap-FFT - FFT direta'); grid on;

nexttile;                                                   % posição 3ª linha, 3ª coluna 
bar([t_eq, t_conv, t_fft, t_ol_conv, t_ol_fft]);                   % barra comparativa dos tempos médios
set(gca,'XTickLabel',{'EqDif','Conv','FFT','OL-Conv','OL-FFT'});   % rótulos das barras
ylabel('Tempo (s)'); title('Comparação de tempo de execução');

% ----------------- Reprodução dos sinais filtrados (Overlap-Add) ---
tocou = false;                                                      % flag para saber se tocou áudio
choice = input('Deseja ouvir o áudio após filtragem por Convolução Overlap-add? (s/n): ','s'); % pergunta ao usuário
if lower(choice) == 's'                                             % se escolher 's' (sim)
    fprintf('Reproduzindo áudio filtrado (Convolução Overlap-Add)...\n'); % informa
    sound(y_ol_conv, fs);                                           % toca y_ol_conv com frequência fs
    pause(length(y_ol_conv)/fs);                                    % pausa até terminar de tocar
    tocou = true;                                                   % marca que tocou
end

choice = input('Deseja ouvir o áudio após filtragem por FFT Overlap-Add? (s/n): ','s');  % pergunta para OL-FFT
if lower(choice) == 's'                                             % se sim
    fprintf('Reproduzindo áudio filtrado (FFT Overlap-Add)...\n');  % informa
    sound(y_ol_fft, fs);                                            % toca y_ol_fft
    pause(length(y_ol_fft)/fs);                                     % pausa até terminar
    tocou = true;                                                   % marca que tocou
end

if ~tocou                                                           % se nenhum áudio foi reproduzido
    fprintf('Ok, nenhum áudio será reproduzido.\n');                % informa que nada será reproduzido
end