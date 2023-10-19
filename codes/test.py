import numpy as np
import matplotlib.pyplot as plt

# Criar um sinal de exemplo
fs = 1000  # Frequência de amostragem (Hz)
t = np.linspace(0, 1, fs, endpoint=False)  # Vetor de tempo de 1 segundo

freq = 5  # Frequência do sinal (Hz)
signal = np.sin(2 * np.pi * freq * t)

# Calcular a FFT do sinal
fft_result = np.fft.fft(signal)

# Calcular as frequências correspondentes à FFT
frequencies = np.fft.fftfreq(len(fft_result), 1 / fs)


# Comprimento do sinal
N = len(fft_result)

# Intervalo de amostragem no tempo
delta_t = 100.0  # Substitua pelo valor apropriado em seu caso

# Inicializa um array vazio para armazenar a DFT inversa
ndiv10 = int(N/10)
idft_result = np.zeros(ndiv10, dtype=complex)

tr = np.linspace(0, 1, int(fs/10), endpoint=False)  # Vetor de tempo de 1 segundo
# Calcula a DFT inversa manualmente e relaciona ao tempo
for n in range(ndiv10):
    t = n*10
    for k in range(N):
        idft_result[n] += fft_result[k] * np.exp(1j * 2 * np.pi * k * t / (N)) / N


N = len(fft_result)
n = len(tr)
ifft_result = np.zeros(n, dtype=complex)
for i in range(n):
    p = i/n
    for k in range(N):
        ifft_result[i] += fft_result[k] * np.exp (1j * 2 * np.pi * k * p) / N

# print(frequencies)
# # Plotar o sinal original e sua derivada
# plt.figure(figsize=(12, 6))
# plt.plot(frequencies)
# plt.title('Frequencies')
# plt.legend()
# plt.show()
t = np.linspace(0, 1, fs, endpoint=False)  # Vetor de tempo de 1 segundo

# Plotar o sinal original e sua derivada
plt.figure(figsize=(12, 6))
plt.plot(t, signal, label="signal")
plt.plot(tr, idft_result, '.', label="fourier")
plt.plot(tr, ifft_result, '.', label="My fourier")
# plt.plot(t, np.fft.ifft(fft_result), ',', label="fourier")
plt.title('Sinal Original')
plt.legend()
plt.show()

