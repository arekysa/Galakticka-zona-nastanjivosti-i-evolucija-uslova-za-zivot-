import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import skew, kurtosis
from scipy.optimize import curve_fit

# globalni parametri za stilizaciju grafika
rc('text', usetex=True)
rc('font', size=12)
rc('axes', titlesize=16)
rc('axes', labelsize=16)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)
rc('legend', fontsize=12)
rc('figure', titlesize=16)
plt.rcParams["font.family"] = "Times New Roman"


# ucitavanje podataka
df = pd.read_csv('mw_analogs.csv')

# prebacivanje u log10 bolje za grafike
df['log_Mstar'] = np.log10(df['StellarMass'])
df['log_SFR'] = np.log10(df['SFR'])
df['log_Zgas'] = np.log10(df['GasMetallicity'])
df['log_Zstar'] = np.log10(df['StarMetallicity'])

# specificna stopa formiranja zvezda
# df['log_sSFR'] = df['log_SFR'] - df['log_Mstar']

print(f"Uspesno ucitano i obradjeno {len(df)} galaksija.")
print(df.head()) # prikaz prvih 5 cisto za proveru


def statisticka_analiza(kolona, promenljiva):
    
    podaci = df[kolona].replace([np.inf, -np.inf], np.nan) # ako negde postoje beskonacnosti pretvara ih u NaN
    podaci = df[kolona].dropna() # cistimo podatke od NaN vrednosti
    
    aritmeticka_sredina = np.mean(podaci)
    standardna_devijacija = np.std(podaci, ddof=1)
    
    # Robusne mere
    medijana = np.median(podaci)
    q75, q25 = np.percentile(podaci, [75, 25])
    sigma_G = q75 - q25 # IQR (interkvartilni opseg)
    
    # oblik raspodele
    sk = skew(podaci)
    kurt = kurtosis(podaci)
    
    print(f'--- {promenljiva} ---')
    print(f'  Srednja vr. (Mean) : {aritmeticka_sredina:.3f}')
    print(f'  Medijana (Median)  : {medijana:.3f}')
    print(f'  Stand. dev. (Std)  : {standardna_devijacija:.3f}')
    print(f'  Sigma_G (IQR)      : {sigma_G:.3f}')
    print(f'  Skewness           : {sk:.3f}')
    print(f'  Kurtosis           : {kurt:.3f}\n')
    
    return aritmeticka_sredina, medijana

print('=== Statisticka analiza uzorka ===\n')
mean_m, med_m = statisticka_analiza('log_Mstar', 'Stellar Mass (log M*)')
mean_sfr, med_sfr = statisticka_analiza('log_SFR', 'Star Formation Rate (log SFR)')
mean_z, med_z = statisticka_analiza('log_Zgas', 'Gas Metallicity (log Z_gas)')

"""
=== Statisticka analiza uzorka ===

--- Stellar Mass (log M*) ---
  Srednja vr. (Mean) : 10.646
  Medijana (Median)  : 10.626
  Stand. dev. (Std)  : 0.167
  Sigma_G (IQR)      : 0.272
  Skewness           : 0.384
  Kurtosis           : -0.946

--- Star Formation Rate (log SFR) ---
  Srednja vr. (Mean) : -0.200
  Medijana (Median)  : 0.131
  Stand. dev. (Std)  : 0.981
  Sigma_G (IQR)      : 1.136
  Skewness           : -1.268
  Kurtosis           : 1.027

--- Gas Metallicity (log Z_gas) ---
  Srednja vr. (Mean) : -1.781
  Medijana (Median)  : -1.760
  Stand. dev. (Std)  : 0.187
  Sigma_G (IQR)      : 0.257
  Skewness           : -0.320
  Kurtosis           : 0.012

"""




fig, axes = plt.subplots(1, 3, figsize=(16, 4))
boje = ['#3B82F6', '#22C55E', '#EF4444']
kolone = ['log_Mstar', 'log_SFR', 'log_Zgas']
naslovi = [r'Zvezdana masa', r'Stopa formiranja zvezda', r'Metalicnost gasa']
x_ose = [r'$\log(M_\ast [M_\odot])$', r'$\log(\rm{SFR})$', r'$\log(Z_{gas})$']

for i in range(3):
    podaci = df[kolone[i]].dropna()
    mean_val = np.mean(podaci)
    med_val = np.median(podaci)
    
    axes[i].hist(podaci, bins=60, color=boje[i], alpha=0.8, edgecolor='white')
    
    axes[i].axvline(mean_val, color='black', ls='--', lw=2, label=f'Mean: {mean_val:.2f}')
    axes[i].axvline(med_val, color='black', ls=':', lw=2, label=f'Median: {med_val:.2f}')
    axes[i].set_title(naslovi[i])
    axes[i].set_xlabel(x_ose[i])
    axes[i].set_ylabel('Broj galaksija')
    axes[i].legend()

plt.tight_layout()
plt.savefig('medijana_i_srv.pdf')
plt.show()





