%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:
%   - Brian Nicholas Alves Rondon Lahr
%   - Fernando Cirilo Zanchetta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Configurações do Estudo
% -------------------------------------------------------------------------
% Quantidade de bits/níveis de quantização de fase do holograma final
phase_quantization_bits_to_test = 1:16;
num_tests = length(phase_quantization_bits_to_test);

% Vetores para guardar os dados e plotar depois
vec_eff_theo = zeros(1, num_tests);
vec_eff_sim  = zeros(1, num_tests);
vec_snr      = zeros(1, num_tests);

% Cria a pasta para salvar as imagens se ela não existir
output_folder = 'hologramas';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Variável para controle do salvamento de imagens (default: true)
saveImages = true;

% -------------------------------------------------------------------------
% Parâmetros Físicos
% -------------------------------------------------------------------------
lambda = 633e-9;      % Comprimento de Onda (Laser Vermelho HeNe 633nm)
pixel_pitch = 8e-6;   % Pixel do Holograma/SLM (8 microns)
z_prop = 0.3;         % Distância de Propagação
k0 = 2 * pi / lambda; % Número de Onda

% Número de pixels em cada direção
Nx = 1080; Ny = 1080;

% Parâmetros do IFTA
k_iter = 25;          % Número máximo de iterações
alpha = 1.1;          % Fator de reforço da imagem
psi = k_iter;         % Psi configurado para o limite superior de c

% -------------------------------------------------------------------------
% Geometria de Fresnel
% -------------------------------------------------------------------------
L_holo = Nx * pixel_pitch;                % Tamanho físico do holograma
pixel_img = (lambda * z_prop) / L_holo;   % Tamanho do pixel no plano de reconstrução
L_img = Nx * pixel_img;                   % Tamanho total da imagem na plano de reconstrução

fprintf('Tamanho do Holograma: %.6f mm\n', L_holo * 1000);
fprintf('Tamanho da Imagem Projetada: %.6f mm (em z = %.6fm)\n', L_img * 1000, z_prop);

% Define o espaço baseado no tamanho do holograma
x = linspace(-L_holo/2, L_holo/2, Nx);
[X, Y] = meshgrid(x, x);

% Fator de Fase Quadrática
Q_phase = exp(1i * k0 / (2 * z_prop) * (X .^ 2 + Y .^ 2));
Q_phase = fftshift(Q_phase);

% -------------------------------------------------------------------------
% Leitura da Imagem Alvo
% -------------------------------------------------------------------------
img = imread('USP-Brasão.png');
img = imcomplement(img);
img = im2double(img);
img = img(:,:,1);
img = imresize(img,[Nx Ny]);
img = img / max(img(:));

% Definição da Região de Interesse (ROI / Omega)
Omega = img > 0.005;

figure;
colormap gray; colorbar;
subplot(1,2,1); imagesc(img); axis image; title('Imagem Alvo');
subplot(1,2,2); imagesc(Omega); axis image; title('Máscara Omega');

% -------------------------------------------------------------------------
% Inicialização
% -------------------------------------------------------------------------
fprintf('\nIniciando Varredura de Bits...\n');

for bits_to_test_counter = 1:num_tests

  current_bits = phase_quantization_bits_to_test(bits_to_test_counter);
  final_levels = 2 ^ current_bits;

  fprintf('\n--- Testando %d Bits (%d Níveis) ---\n', current_bits, final_levels);

  fprintf('Iniciando IFTA com Fresnel...\n');

  % Chute de uma fase inicial no holograma
  g_hologram = exp(1i * 2 * pi * rand(Nx, Ny));

  % Inicializa uma dummy variable e uma dummy matrix para armazenar futuramente a melhor SNR e o melhor holograma
  bestSNR = -inf;
  bestH = zeros(Nx, Ny);

  % -------------------------------------------------------------------------
  % Loop IFTA Híbrido
  % -------------------------------------------------------------------------


  for c = 1 : psi

      % Reduz os níveis de quantização gradualmente
      levels_c = round(final_levels * (psi - c + 1) / psi);
      levels_c = max(levels_c, 2);

      fprintf('  Iteração %d/%d: Níveis de Fase = %d', c, psi, levels_c);

      for iter = 1:k_iter % Sub-loop para estabilizar cada nível

          % --- Propagação FRENTE (Holograma -> Plano de Reconstrução) ---

          % Quantização no plano do holograma
          phase_h = angle(g_hologram);
          phase_q = round((phase_h + pi) / (2 * pi) * (levels_c - 1));
          phase_q = phase_q / (levels_c - 1) * 2 * pi - pi;

          % Amplitude Uniforme
          U_h_quantized = 1.0 .* exp(1i * phase_q);

          U_reconst = fftshift(fft2(fftshift(U_h_quantized .* Q_phase)));

          % --- Restrições no Plano da Imagem (Plano de Reconstrução) ---

          % Amplitude
          Amp_rec = abs(U_reconst);
          Phase_rec = angle(U_reconst);

          % Novo Campo Alvo
          g_prime = zeros(Nx, Ny);

          % DENTRO de Omega (Existe Imagem):
          % Atualização da Amplitude de Entrada
          g_prime(Omega) = alpha * img(Omega) .* exp(1i * Phase_rec(Omega));

          % FORA de Omega (Fundo Preto):
          if c < psi * 0.8 % Nas primeiras iterações, força zero fora para concentrar energia
               g_prime(~Omega) = 0;
          else
               % Nas finais, deixa o ruído existir fora para limpar dentro
               g_prime(~Omega) = U_reconst(~Omega);
          end

          % --- Propagação TRÁS (Plano de Reconstrução -> Holograma) ---
          U_back = ifftshift(ifft2(ifftshift(g_prime)));
          g_hologram = U_back .* conj(Q_phase);
      end

      % -------------------------------
      % Avaliação da SNR
      % -------------------------------
      I_rec = abs(U_reconst) .^ 2;
      scale = max(I_rec(Omega));
      I_rec = I_rec / scale;
      signal = sum(I_rec(Omega) .^ 2);
      noise = sum(I_rec(~Omega) .^ 2);

      mse = sum((I_rec(Omega) - img(Omega)) .^ 2 ) / sum(Omega(:));
      SNR = 10 * log10(1 / mse);

      % -------------------------------
      % Critério de Parada / Armazenamento
      % -------------------------------
      if SNR > bestSNR
          bestSNR = SNR;
          bestH = phase_q;
      elseif SNR < 0.95 * bestSNR
          fprintf(' - SNR: %.6f dB\n\n', SNR);
          fprintf('  Convergência Atingida.\n\n');
          break;
      end

      fprintf(' - SNR: %.6f dB\n\n', SNR);
  end


  % -------------------------------------------------------------------------
  % Resultado Final
  % -------------------------------------------------------------------------

  % Usa Melhor Holograma na Imagem Final
  U_final_h = exp(1i * bestH);
  U_final_img = fftshift(fft2(fftshift(U_final_h .* Q_phase)));
  I_final = abs(U_final_img) .^ 2;

  % Normalização
  I_show = I_final / max(I_final(:));

  if saveImages
    fig = figure('Visible', 'off', 'Position', [0 0 1200 500]);

    axis_img = linspace(-L_img/2, L_img/2, Nx) * 1000;

    subplot(1,2,1);
    imagesc(bestH); axis image; colormap gray;
    title(sprintf('Holograma de Fase (%d Bits)', current_bits));
    xlabel('Pixel x'); ylabel('Pixel y');

    subplot(1,2,2);
    imagesc(axis_img, axis_img, I_show);
    axis image; colormap gray;
    colorbar;
    title(sprintf('Reconstrução (z = %.2fm)(SNR: %.6fdB)', z_prop, bestSNR));
    xlabel('mm'); ylabel('mm');

    filename = sprintf('%s/Resultado_%02d_Bits.png', output_folder, current_bits);
    print(fig, filename, '-dpng', '-r600');

    fprintf('   -> Figura salva: %s\n', filename);
    close(fig); % Fecha a figura da memória
  else
    figure;
    axis_img = linspace(-L_img/2, L_img/2, Nx) * 1000;

    subplot(1,2,1);
    imagesc(bestH); axis image; colormap gray;
    title('Holograma de Fase Final');

    subplot(1,2,2);
    imagesc(axis_img, axis_img, I_show);
    axis image; colormap gray;
    title(sprintf('Reconstrução IFTA + Fresnel (z = %.2fm)', z_prop));
    xlabel('mm'); ylabel('mm');
    colorbar;
  end

  % -------------------------------------------------------------------------
  % Eficiências
  % -------------------------------------------------------------------------

  % Eficiência Teórica (Baseada nos Níveis de Fase)
  eta_teorica = (sinc(1 / final_levels)) ^ 2;

  % Eficiência Simulada (Energia no Sinal / Energia Total)
  Energy_Total = sum(I_final(:));
  Energy_Signal = sum(I_final(Omega));
  eta_simulada = Energy_Signal / Energy_Total;

  fprintf('\n--- Resultados de Eficiência ---\n');
  fprintf('Eficiência Teórica (%d níveis): %.6f%%\n', final_levels, eta_teorica * 100);
  fprintf('Eficiência Simulada: %.6f%%\n', eta_simulada * 100);
  fprintf('SNR Final: %.6f dB\n', bestSNR);

  % Armazena nos Vetores de Análise
  vec_eff_theo(bits_to_test_counter) = eta_teorica;
  vec_eff_sim(bits_to_test_counter)  = eta_simulada;
  vec_snr(bits_to_test_counter)      = bestSNR;
end

% Geração dos Gráficos de Eficiência
subplot(1, 2, 1);
plot(phase_quantization_bits_to_test, vec_eff_theo * 100, '-o', 'LineWidth', 2, 'DisplayName', 'Teórica');
hold on;
plot(phase_quantization_bits_to_test, vec_eff_sim * 100, '-s', 'LineWidth', 2, 'DisplayName', 'Simulada');
grid on;
xlabel('Número de Bits (Quantização)');
ylabel('Eficiência (%)');
title('Eficiência de Difração');
legend('Location', 'SouthEast');
ylim([0 105]);

% Gráfico 2: SNR
subplot(1, 2, 2);
plot(phase_quantization_bits_to_test, vec_snr, '-d', 'LineWidth', 2, 'Color', '#77AC30');
grid on;
xlabel('Número de Bits');
ylabel('SNR (dB)');
title('Qualidade da Imagem (SNR)');
