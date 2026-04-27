# kod koji izvlaci analoge MW iz TNG100-1 iz snap=99 (z=0) 

import illustris_python as il
import pandas as pd
import numpy as np

basePath = "/home/arekysa/Documents/Praksa/AOB_Ana_M_24_25/TNG100-1"
snap = 99
h = 0.6774

# ucitavanje specificnik polja iz kataloga
polja = ['SubhaloMassType', 'SubhaloSFR', 'SubhaloFlag']
subhalos = il.groupcat.loadSubhalos(basePath, snap, fields=polja)

# ekstrakcija mase zvezda
stellar_mass = subhalos['SubhaloMassType'][:, 4] * 1e10 / h
sfr = subhalos['SubhaloSFR']
flag = subhalos['SubhaloFlag']

# definisanje granica za selekciju i filtriranje
mass_min = 10**10.4
mass_max = 10**11.0  
sfr_min = 0.5
sfr_max = 3.0 # prvobitno stvaljeno 5 i dobijeno je 968 galaksija


mw_condition = (stellar_mass >= mass_min) & \
               (stellar_mass <= mass_max) & \
               (sfr >= sfr_min) & \
               (sfr <= sfr_max) & \
               (flag == 1)

mw_ids = np.where(mw_condition)[0]

# data frame za konverziju u csv
mw_analogs = pd.DataFrame({
    'SubhaloID': mw_ids,
    'StellarMass': stellar_mass[mw_ids],
    'SFR': sfr[mw_ids]
})

print(f"Pronađeno je {len(mw_ids)} MW analoga.")
output_path = "mw_analogs_ids_2.csv"
mw_analogs.to_csv(output_path, index=False)

