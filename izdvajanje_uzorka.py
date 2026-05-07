# kod koji izvlaci analoge MW iz TNG100-1 iz snap=99 (z=0) 
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
import illustris_python as il
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np

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


basePath = "/home/arekysa/Documents/Praksa/AOB_Ana_M_24_25/TNG100-1"
snap = 99
h = 0.6774

# DODATO: 'SubhaloStarMetallicity'
polja = ['SubhaloMassType', 'SubhaloSFR', 'SubhaloGasMetallicity', 'SubhaloStarMetallicity', 'SubhaloFlag']
subhalos = il.groupcat.loadSubhalos(basePath, snap, fields=polja)

# ekstrakcija masa
stellar_mass = subhalos['SubhaloMassType'][:, 4] * 1e10 / h
gas_mass = subhalos['SubhaloMassType'][:, 0] * 1e10 / h  
z_gas = subhalos['SubhaloGasMetallicity']
z_star = subhalos['SubhaloStarMetallicity']             
sfr = subhalos['SubhaloSFR']
flag = subhalos['SubhaloFlag']

# definisanje granica za selekciju i filtriranje
mass_min = 10**10.4
mass_max = 10**11.0  

mw_condition = (stellar_mass >= mass_min) & \
               (stellar_mass <= mass_max) & \
               (sfr > 0) & \
               (z_gas > 0) & \
               (flag == 1)

mw_ids = np.where(mw_condition)[0]

# data frame za konverziju u csv
mw_analogs = pd.DataFrame({
    'SubhaloID': mw_ids,
    'StellarMass': stellar_mass[mw_ids],
    'GasMass': gas_mass[mw_ids],          
    'SFR': sfr[mw_ids],
    'GasMetallicity': z_gas[mw_ids],
    'StarMetallicity': z_star[mw_ids]  
})

print(f'Pronađeno je {len(mw_ids)} MW analoga.') # Pronađeno je 1703 MW analoga.
output_path = 'mw_analogs.csv'
mw_analogs.to_csv(output_path, index=False)


#################### grafik globalne price i iskljucivo nase populacije #################################

# uzimamo polja koja nam trebaju za grafik
fields = ['SubhaloFlag', 'SubhaloMassType', 'SubhaloSFR', 'SubhaloStarMetallicity']
subs = il.groupcat.loadSubhalos(basePath, 99, fields=fields)
subhalos = pd.DataFrame()
subhalos['SubhaloFlag'] = subs['SubhaloFlag']
subhalos['log_Mstar'] = np.log10(subs['SubhaloMassType'][:, 4] * 1e10 / h)
subhalos['SFR'] = subs['SubhaloSFR']
subhalos['log_Zstar'] = np.log10(subs['SubhaloStarMetallicity'])

# filtriranje
mask_broad = (
    (subhalos['log_Mstar'] >= 8) & 
    (subhalos['SFR'] > 0) & 
    np.isfinite(subhalos['log_Zstar']) & 
    (subhalos['SubhaloFlag'] == 1)
)
df_broad = subhalos[mask_broad]

plt.figure(figsize=(10, 6))
plt.hist2d(df_broad['log_Mstar'], df_broad['log_Zstar'], bins=150, cmap='viridis', cmin=1, norm=LogNorm(vmin=1))
plt.axvline(10.4, color='red', ls='--', label='Izabarani opseg (10.4 do 11.0)')
plt.axvline(11.0, color='red', ls='--')
plt.axvspan(10.4, 11.0, color='red', alpha=0.05)
plt.xlim(8, 12)
plt.xlabel(r'$\log(M_\ast  [M_\odot])$')
plt.ylabel(r'$\log(Z_\ast)$')
plt.title('Relacija Masa-Metalicnost (TNG100) — Pozicija MW analoga')
plt.colorbar(label='Broj galaksija', pad=0.01)
plt.legend(loc='lower left', framealpha=0.8)
plt.tight_layout()
plt.tick_params(direction='in', top=True, right=True)
plt.savefig('relacija_masa_metalicnost.pdf')
plt.show()

mask_zoom = (df_broad['log_Mstar'] >= 10.4) & (df_broad['log_Mstar'] <= 11.0)
df_zoom = df_broad[mask_zoom]

def linearni_fit(x, k, n):
    return k * x + n

mask_fit = np.isfinite(df_zoom['log_Mstar']) & np.isfinite(df_zoom['log_Zstar'])
x_podaci = df_zoom.loc[mask_fit, 'log_Mstar'].values
y_podaci = df_zoom.loc[mask_fit, 'log_Zstar'].values

popt, pcov = curve_fit(linearni_fit, x_podaci, y_podaci)
k_fit, n_fit = popt

print(f'=== REZULTATI FITA ===')
print(f'Nagib (k)   = {k_fit:.3f}')
print(f'Odsecak (n) = {n_fit:.3f}\n')

"""
=== REZULTATI FITA ===
Nagib (k)   = 0.041
Odsecak (n) = -2.084

"""

plt.figure(figsize=(8, 6))
plt.hist2d(df_zoom['log_Mstar'], df_zoom['log_Zstar'], bins=50, cmap='viridis', cmin=1, norm=LogNorm(vmin=1))

x_linija = np.linspace(10.4, 11.0, 100)
y_linija = linearni_fit(x_linija, k_fit, n_fit)
plt.plot(x_linija, y_linija, color='red', lw=3, label=f'Linearni fit\n$y = {k_fit:.2f}x {n_fit:+.2f}$')
plt.xlabel(r'$\log(M_\ast [M_\odot])$')
plt.ylabel(r'$\log(Z_\ast)$')
plt.title('Relacija Masa-Metalicnost za MW analoge sa linearnim fitom')
plt.colorbar(label='Broj galaksija', pad=0.01)
plt.legend(loc='lower right', framealpha=0.8)
plt.tight_layout()
plt.savefig('zumirano_masa_metalicnost.pdf')
plt.show()
