import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 10

def norm_min_max(x):
    min_val = np.min(x)
    max_val = np.max(x)
    x_norm = (x - min_val) / (max_val - min_val)
    return x_norm

def norm_PSD(fp, x):
    S = np.fft.fft(x)
    Pa = np.sqrt(S*np.conj(S))
    Pa = Pa.real / len(x)
    f = np.fft.fftfreq(len(x), d=1/fp)
    Pa = Pa[np.logical_and(f >= 0, f <= 12000)]
    f = f[np.logical_and(f >= 0, f <= 12000)]
    Pa_norm = norm_min_max(Pa)
    return f, Pa_norm

def norm_welch_PSD(fp, x):
    (f, Pa) = scipy.signal.welch(x, fp, window='hann', scaling='spectrum', average='median', nperseg=2500)
    Pa = Pa[np.logical_and(f >= 0, f <= 12000)]
    f = f[np.logical_and(f >= 0, f <= 12000)]
    Pa_norm = norm_min_max(Pa)
    return f, Pa_norm

def real_PSD(fp, x):
    S = np.fft.rfft(x)
    Pa = np.sqrt(S*np.conj(S))
    Pa = Pa.real / len(x)
    f = np.fft.rfftfreq(len(x), d=1/fp)
    Pa = Pa[np.logical_and(f >= 0, f <= 12000)]
    f = f[np.logical_and(f >= 0, f <= 12000)]
    return f, Pa
  
# --- Basic gearbox data

z1=30
z2=47
ns=38 * 60 / 2.5  # prędkość obrotowa silnika
n1=2279           # prędkość obrotowa zębnika
n2=n1*z1/z2       # prędkość obrotowa koła
fs=ns/60
f1=n1/60
f2=n2/60
fz=f1*z1
fk=fz/np.lcm(z1, z2)

# --- Vibration measurement data preprocessing

dataframe = pd.read_csv('DRG.csv', header=None)
DRG = dataframe.values
fp = 25000

i_r, i_c = DRG.shape
a_meas = np.reshape(DRG, i_r * i_c)
t_meas = np.linspace(0, i_r, i_r*fp)

(f, Pa) = norm_PSD(fp, a_meas)
(f_w, Pa_w) = norm_welch_PSD(fp, a_meas)
(f_r, Pa_r) = real_PSD(fp, a_meas)

am_peaks, _ = scipy.signal.find_peaks(Pa, height=0.1)

# --- Transmission error preprocessing

TEdataframe = pd.read_csv('TE_min2.csv', header=None)
TE = TEdataframe.values

t = TE[0, :]
x = TE[1, :] * 10e6 # na mikrometry
v = TE[2, :]
a = TE[3, :]
dt = np.diff(t)
fp_te = 1 / dt[0]

(fx, Pxx) = norm_PSD(fp_te, x)
(fa, Paa) = norm_PSD(fp_te, a)

(fx_w, Pxx_w) = norm_welch_PSD(fp_te, x)
(fa_w, Paa_w) = norm_welch_PSD(fp_te, a)

(fx_r, Pxx_r) = real_PSD(fp_te, x)
(fa_r, Paa_r) = real_PSD(fp_te, a)

x_peaks, _ = scipy.signal.find_peaks(Pxx, height=0.1)
a_peaks, _ = scipy.signal.find_peaks(Paa, height=0.1)

# --- Interpolation

t_int = np.linspace(0,1, 25000)
x_int = np.interp(t_int, t, x)
a_int = np.interp(t_int, t, a)

# --- Measurement data visualisation

fig, axs = plt.subplots(ncols=1, nrows=2)
axs[0].plot(t_meas[t_meas <= 0.025], a_meas[t_meas <= 0.025])
axs[0].set_xlabel('Time, s')
axs[0].set_ylabel('a, $\mathdefault{m/s^2}$')
axs[0].set_title('a)')
axs[0].grid(visible=True)
#axs[0].autoscale(enable=True, axis='x', tight=True)
axs[0].autoscale(enable=True, axis='y', tight=True)
axs[0].set_xlim(left=0, right=0.02)

axs[1].plot(f_r, Pa_r)
axs[1].set_xlabel('Frequency, Hz')
axs[1].set_ylabel('A(a) $\mathdefault{m/s^2}$')
axs[1].set_title('b)')
axs[1].grid(visible=True)
axs[1].autoscale(enable=True, axis='x', tight=True)
axs[1].autoscale(enable=True, axis='y', tight=True)

fig.tight_layout()
plt.show()

# --- Visualisation of transmission error excitation

fig, axs = plt.subplots(3, 1)

ax = axs[0]
ax.plot(t[np.logical_and(t >=0.002, t <= 0.006)], x[np.logical_and(t >=0.002, t <= 0.006)])
ax.set_xlabel('Time, s')
ax.set_ylabel('$\mathdefault{x_{te}}$, $\mathdefault{\mu m}$')
ax.set_title('a)')
ax.set_xlim(left=0.002, right=0.006)
ax.grid(visible=True)
ax.autoscale(enable=True, axis='x', tight=True)
ax.autoscale(enable=True, axis='y', tight=True)

ax = axs[1]
ax.plot(fx, Pxx)
ax.set_xlabel('Frequency, Hz')
ax.set_ylabel('$\mathdefault{A_{norm}(x_{te})}$')
ax.set_title('b)')
ax.grid(visible=True)
ax.autoscale(enable=True, axis='x', tight=True)
ax.autoscale(enable=True, axis='y', tight=True)

ax = axs[2]
ax.plot(fa, Paa)
ax.set_xlabel('Frequency, Hz')
ax.set_ylabel('$\mathdefault{A_{norm}(a_{te})}$')
ax.set_title('c)')
ax.grid(visible=True)
ax.autoscale(enable=True, axis='x', tight=True)
ax.autoscale(enable=True, axis='y', tight=True)

fig.tight_layout()
plt.show()

# --------- Correlation with measurement data

fax, Pax = scipy.signal.csd(norm_min_max(a_meas), norm_min_max(x_int), fs=25000,  nperseg=1250)
Pax = norm_min_max(np.abs(Pax))
ax_peaks, _ = scipy.signal.find_peaks(Pax)
fr_ax = fax[ax_peaks]

faaa, Paaa = scipy.signal.csd(norm_min_max(a_meas), norm_min_max(a_int), fs=25000,  nperseg=1250)
Paaa = norm_min_max(np.abs(Paaa))
aaa_peaks, _ = scipy.signal.find_peaks(Paaa)
fr_aaa = fax[aaa_peaks]

fig, axs = plt.subplots(2, 1)

ax = axs[0]
ax.plot(fax, Pax)
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('$\mathdefault{CSD(a, x_{te})}$')
ax.grid(visible=True)
ax.autoscale(enable=True, axis='x', tight=True)
ax.set_xlim(left=0, right=12000)
ax.set_title('a)')
x_labels = []
x_tic = np.zeros(9)
for i in range(9):
    x_tic[i] = np.floor((i + 1) * fz)
    x_labels.append(str(i+1) + '$\mathdefault{f_m=}$\n' + '%g'%(x_tic[i]))

ax.set_xticks(ticks=x_tic, labels=x_labels)
ax.autoscale(enable=True, axis='x', tight=True)
ax.autoscale(enable=True, axis='y', tight=True)

ax = axs[1]
ax.plot(faaa, Paaa)
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('$\mathdefault{CSD(a, a_{te})}$')
ax.grid(visible=True)
ax.autoscale(enable=True, axis='x', tight=True)
ax.set_xlim(left=0, right=12000)
ax.set_title('b)')
ax.set_xticks(ticks=x_tic, labels=x_labels)
ax.autoscale(enable=True, axis='x', tight=True)
ax.autoscale(enable=True, axis='y', tight=True)

fig.tight_layout()
plt.show()
