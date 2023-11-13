import numpy as np
from scipy.fft import nufft

# Exemplo de amostras não uniformemente espaçadas
x_original = np.array([1, 2, 3, 4])
t = np.array([0, 1, 2, 3])

# Frequências associadas aos índices de saída
frequencies = np.fft.fftfreq(len(x_original))

# Calcula a NUDFT
X_nufft = nufft(x_original, t, frequencies)

# Reconstroi o sinal original usando a NUDFT inversa
N = len(x_original)
x_reconstructed = np.zeros(N, dtype=np.complex128)

for n in range(N):
    x_reconstructed[n] = np.sum(X_nufft * np.exp(1j * 2 * np.pi * frequencies * t[n]))

# Imprime os resultados
print("Sinal Original:", x_original)
print("Sinal Reconstruído:", np.real(x_reconstructed))  # Mostra a parte real, pois o sinal é real
