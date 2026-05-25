#!/usr/bin/env python
"""Tárgy: Együttműködési lehetőség — AI-alapú gyógyszer-újrahasznosítási platform"""

KEDVEZMÉNYEZETT = "Semmelweis Egyetem, Farmakológiai és Farmakoterápiás Intézet"
TÁRGY = "Együttműködési javaslat: Mesterséges intelligencia alapú gyógyszer-újrahasznosítási platform validációja"

LEVÉL = f"""
Tisztelt Intézetvezető Úr/Asszony!

Egy egyetemi hallgatói projekt keretében kifejlesztettünk egy nyílt forráskódú, 
mesterséges intelligenciával támogatott gyógyszer-újrahasznosítási (drug repurposing) 
platformot, amely képes tetszőleges fehérje célpont ellen automatikusan szűrni 
FDA által jóváhagyott gyógyszereket és azok optimalizált variánsait.

A platform jelenlegi eredményei:

1. Nilotinib_Var_17 – egy optimalizált Nilotinib variáns – az alábbi 
   fehérje célpontokon mutat kiemelkedő in silico kötődést:

   • SARS-CoV-2 Mpro (COVID-19)         -8.85 kcal/mol  (PDB: 6LU7)
   • HIV-1 Proteáz                      -10.09 kcal/mol  (PDB: 1HPV)
   • EGFR Kináz (daganat)               -10.94 kcal/mol  (PDB: 1M17)
   • Bcl-2 (apoptózis)                  -8.06 kcal/mol   (PDB: 4LVT)
   • COX-2 (gyulladás)                  erős kötődés     (PDB: 5KIR)

2. A platform automatikusan:
   • Letölti és előkészíti a fehérje szerkezetet (RCSB PDB)
   • Detektálja a kötőzsebet
   • AI predikcióval előszűr 2000+ FDA gyógyszert
   • AutoDock Vina molekuláris dokkolást végez
   • Rangsorolt riportot készít

3. A gépi tanulási modell R²=0.63 pontossággal prediktálja a kötési affinitást
   5 különböző fehérjén végzett validáció alapján.

A projekt nyílt forráskódú (MIT licenc), elérhető GitHub-on:
https://github.com/horvatjanos/universal-drug-repurposing

Keresünk akadémiai partnert a prediktált gyógyszerjelöltek in vitro validációjához. 
Az együttműködés keretében az alábbiakat tudjuk biztosítani:

• A platform ingyenes és korlátlan használata
• Testreszabás az Intézet specifikus kutatási céljaihoz
• Technikai támogatás és továbbfejlesztés
• Közös publikáció lehetősége

Amennyiben érdeklődnek az együttműködés iránt, kérjük jelezzék a fenti GitHub 
oldalon vagy az alábbi elérhetőségen.

Tisztelettel,

Horváth János
E-mail: [email]
GitHub: https://github.com/horvatjanos
Projekt: https://github.com/horvatjanos/universal-drug-repurposing

---

Melléklet: A Nilotinib_Var_17 molekula SMILES kódja és prediktált tulajdonságai
"""

print(LEVÉL)
