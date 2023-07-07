import numpy as np
from scipy.signal import butter, filtfilt

def emg_labiocom(dat, freq, graficos, windowlength, overlap):
    # emg_analisys Realiza a análise de EMG de 1 músculo Author: Prof.
    # Dr. Paulo Roberto Pereira Santiago - EEFERP/USP 04/11/2013 %%
        if len(dat) == 1:
        freq = 2000
        windowlength = 3
        overlap = 1.5
        graficos = 1
    elif len(dat) == 2:
        graficos = 1
        windowlength = 0.125
        overlap = 0.0625
    elif len(dat) == 3:
        windowlength = 0.125
        overlap = windowlength/2
    elif len(dat) == 4:
        overlap = 0.0625
    
    windowlength = windowlength * freq
    overlap = overlap * freq

    datemgvoltraw = dat
    datemgvoltfilt = butter_bandpass_filter(datemgvoltraw, 15, 500, freq)

    datemg = datemgvoltfilt[:, 1] * 1e6  # Convertendo de volts para microvolts

def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

import numpy as np
import matplotlib.pyplot as plt

# tempo
nl1, nc1 = datemg.shape
temp1 = np.arange(nl1) / freq

# GRÁFICO 1 - RAW EMG
plt.figure(1)
plt.plot(temp1, dat)
plt.ylabel('Amplitude (V)')
plt.title('RAW EMG')
plt.axis('tight')
plt.grid(True)
plt.show()

print('Deseja selecionar um intervalo de tempo para as análises?')
print('Para "NÃO", digite 0')
print('Para "SIM", digite o tempo inicial e final em segundos, ex. [10, 20]')
print('para uma análise do 10º ao 20º segundo')

tanalise = input('Digite sua opção: ')

if len(tanalise) == 1 and tanalise[0] == 0:
    tini = 0
    tfim = nl1 - 1
elif len(tanalise) == 2 and tanalise[0] == 0:
    tini = 0
    tfim = int(tanalise[1] * freq)
elif len(tanalise) == 2 and tanalise[0] != 0:
    tini = int(tanalise[0] * freq)
    tfim = int(tanalise[1] * freq)

signal = datemg[tini:tfim, :]

# tempo
nl, nc = signal.shape
temp = np.arange(nl) / freq
# GRAFICO 1 - RAW EMG
if graficos == 1:
    plt.figure(1)
    plt.clf()

    plt.plot(temp, signal)
    plt.ylabel('Amplitude (µV)')
    plt.title('EMG SELECT AND FILTER')
    plt.grid(True)
    plt.axis('tight')
    plt.show()

# SIGNAL RECTIFIED
signalD = signal - np.mean(signal)  # detrend
signalDabs = np.abs(signalD)

# GRAFICO 2 - FULL-WAVE RECTIFIED EMG
if graficos == 1:
    plt.figure(2)
    plt.plot(temp, signalDabs)
    plt.ylabel('Amplitude (µV)')
    plt.title('FULL-WAVE EMG')
    plt.axis('tight')
    plt.grid(True)
    plt.show()

# LINEAR ENVELOPE
signalenv = np.abs(butter_lowpass_filter(signalDabs, 50, freq))
signal_integ = np.trapz(signalenv, dx=1/freq)

# GRAFICO 3 - LINEAR ENVELOPE
if graficos == 1:
    plt.figure(2)
    plt.plot(temp, signalenv, 'r', linewidth=2)
    plt.ylabel('Amplitude (µV)')
    plt.title('FULL-WAVE & LINEAR ENVELOPE = ' + str(np.round(signal_integ, 1)) + ' µV.s')
    plt.axis('tight')
    plt.grid(True)
    plt.show()

# RMS do sinal - window length e overlap
signalrms2 = rms2(signal, windowlength, overlap, 'NO')

# Tempo RMS
length_rms2 = len(signalrms2)
temp_rms = np.linspace(0, temp[-1], length_rms2)

# RMS de todo o sinal
signalrms = np.sqrt(np.mean(signal ** 2))

# GRAFICO 4 - RMS do EMG
if graficos == 1:
    plt.figure(3)
    plt.plot(temp_rms, signalrms2)
    plt.plot([temp_rms[0], temp_rms[-1]], [signalrms, signalrms], linewidth=2, linestyle='--', color='r')
    plt.ylabel('RMS (µV)')
    plt.title('EMG RMS = ' + str(np.round(signalrms, 1)) + ' µV')
    plt.xlabel('Tempo (s)')
    plt.axis('tight')
    plt.grid(True)
    plt.show()

# Analise no domínio da frequência
MPF_signal, PEAK_signal, F50_signal, F95_signal, F_signal, P_signal = psd2(signal, freq)

from scipy.signal import spectrogram, butter, filtfilt
from scipy.ndimage import gaussian_filter1d

# GRAFICO 5 - PSD
if graficos == 1:
    plt.figure(4)
    plt.plot(F_signal, P_signal)
    plt.ylabel('Power (dB)')
    plt.title('EMG - PSD')
    plt.xlabel('Frequência (Hz)')
    plt.axis('tight')
    plt.grid(True)
    plt.show()

# Frequencia mediana no dominio no temporal
Fs = freq
t = temp
nfft = 2 ** int(np.ceil(np.log2(windowlength)))
window = np.hanning(int(round(windowlength)))
overlap = int(round(overlap))

signal_S, signal_F, signal_T, signal_P = spectrogram(signal, window, overlap, nfft, Fs)

medianfreqs_signal = np.ones(signal_P.shape[1])

for nn in range(signal_P.shape[1]):
    signal_P_normcumsumpsd = np.cumsum(signal_P[:, nn]) / np.sum(signal_P[:, nn])
    Ind1 = np.where(signal_P_normcumsumpsd <= 0.5)[0][-1]
    sizeInd1 = len(Ind1)

    if sizeInd1 == 0:
        medianfreqs_signal[nn] = 3
    elif sizeInd1 != 0:
        medianfreqs_signal[nn] = signal_F[Ind1]

medianfreqs_signal1 = gaussian_filter1d(medianfreqs_signal, 0.5)
medianfreqs_signalf = medianfreqs_signal1
freqmediani = medianfreqs_signal

tam_mf = len(medianfreqs_signalf)
t_mft = np.linspace(0, temp[-1], tam_mf)  # vetor tempo freq median t

# GRAFICOS 6 - Frequencia Mediana na série temporal
if graficos == 1:
    plt.figure(6)
    plt.plot(t_mft, medianfreqs_signal)
    plt.plot(t_mft, medianfreqs_signalf, color='r', linewidth=1)
    plt.xlabel('Tempo (s)')
    plt.ylabel('Hz')
    plt.title('Frequência Mediana')
    plt.axis('tight')
    plt.grid(True)
    plt.show()

integral = signal_integ
rmsmed = signalrms
freqmedian = F50_signal
psdpeak = PEAK_signal
linenveli = signalenv
rmsi = signalrms2

# DECLARATIONS AND INITIALIZATIONS
# Calculates windowed (over- and non-overlapping) RMS of a signal using the specified windowlength
# y = rms(signal, windowlength, overlap, zeropad)
# signal is a 1-D vector
# windowlength is an integer length of the RMS window in samples
# overlap is the number of samples to overlap adjacent windows (enter 0 to use non-overlapping windows)
# zeropad is a flag for zero padding the end of your data...(0 for NO, 1 for YES)
# ex. y=rms(mysignal, 30, 10, 1).  Calculate RMS with window of length 30 samples, overlapped by 10 samples each, and zeropad the last window if necessary
# ex. y=rms(mysignal, 30, 0, 0).  Calculate RMS with window of length 30 samples, no overlapping samples, and do not zeropad the last window
#
# Author: A. Bolu Ajiboye

delta = windowlength - overlap

# CALCULATE RMS

indices = list(range(0, len(signal), delta))
# Zeropad signal
if len(signal) - indices[-1] + 1 < windowlength:
    if zeropad:
        signal = np.pad(signal, (0, indices[-1] + windowlength - len(signal)), 'constant')
    else:
        indices = indices[:next((i for i, val in enumerate(indices) if val + windowlength - 1 <= len(signal)), None)]

y = np.zeros(len(indices))
# Square the samples
signal = np.square(signal)

index = 0
for i in indices:
    index += 1
    # Average and take the square root of each window
    y[index - 1] = np.sqrt(np.mean(signal[i:i + windowlength]))

import numpy as np
from scipy.signal import welch
from scipy.integrate import trapz

def psd2(x, fs, nfft=1024, window=256, noverlap=None, dflag='mean'):
    if len(x) < 1000:
        nfft = 512
        window = 256
    if noverlap is None:
        noverlap = window // 2

    # Power Spectral Density
    f, p = welch(x, fs, window='hanning', nperseg=window, noverlap=noverlap, nfft=nfft)

    # Mean Power Frequency (MPF)
    mpf = trapz(f * p, f) / trapz(p, f)

    # Peak Frequency
    peak = f[np.argmax(p)]

    # 50% and 95% of PSD
    area = np.cumtrapz(p, f)
    f50_idx = np.argmax(area >= 0.5 * area[-1])
    f50 = f[f50_idx] if f50_idx < len(f) else 0
    f95_idx = np.argmax(area >= 0.95 * area[-1])
    f95 = f[f95_idx] if f95_idx < len(f) else 0

    return mpf, peak, f50, f95, f, p

def agradecer_professor_black():
    mensagem = """
    Querido Professor Black,

    Gostaríamos de expressar nossa sincera gratidão pelas suas aulas de biomecânica.
    Você é tão apaixonado por esse assunto que transforma até as leis da física em uma aventura emocionante!

    Suas explicações cheias de entusiasmo e exemplos divertidos tornaram o aprendizado da biomecânica uma jornada inesquecível.
    Nunca esqueceremos das suas demonstrações ousadas e experimentos loucos!

    Agradecemos por abrir nossas mentes para o fascinante mundo da biomecânica e por nos inspirar a explorar além dos limites convencionais.
    Suas aulas foram verdadeiramente transformadoras e nos ajudaram a enxergar a biomecânica de uma forma totalmente nova.

    E, Professor Black, temos uma confissão a fazer... você nos descobriu! Sim, recorremos à ajuda de um assistente virtual chamado ChatGPT para criar essa mensagem criativa.
    Mas não se preocupe, ainda apreciamos e valorizamos seu estilo único e apaixonado de ensino, mesmo que tenhamos um pequeno "aliado" virtual!

    Obrigado, Professor Black, por sua dedicação, excentricidade e por nos ensinar que a biomecânica é mais do que apenas fórmulas e equações, é uma arte que desafia a gravidade!

    Com gratidão (e um pouco de ajuda da inteligência artificial),
    Seus alunos apaixonados pela biomecânica
    """

    print(mensagem)

# Chamar a função para agradecer ao Professor Black
agradecer_professor_black()

