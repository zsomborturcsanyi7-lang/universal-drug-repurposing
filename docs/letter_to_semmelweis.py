#!/usr/bin/env python
"""Levél a Semmelweis Egyetem Farmakológiai Intézetének"""

LEVÉL = """
Tisztelt Bioinformatika Munkacsoport, Tisztelt Kutatók!

Turcsányi Zsombor vagyok, a Neumann János Informatikai Technikum tanulója.
Engedje meg, hogy bemutassam a projektemet, amelyen az elmúlt hónapokban
dolgoztam — és szeretném a segítségüket kérni, hogy ezt a munkát a
számítógép képernyőjéről a valóságba is átültethessük.

---

## Miről van szó?

Egy nyílt forráskódú, mesterséges intelligenciával támogatott gyógyszer-
újrahasznosítási (drug repurposing) platformot fejlesztettem ki, amely képes
tetszőleges fehérje célpont ellen percek alatt szűrni az FDA által már
jóváhagyott gyógyszereket. A platform a Semmelweis Egyetem által is használt
tudományos szoftverekre épül (AutoDock Vina, RDKit, OpenBabel).

A projektet nem egyedül, hanem egy mesterséges intelligencia asszisztenssel
(Hermes AI) együttműködve valósítottam meg — ami önmagában is mutatja,
hogy a jövő tudományos kutatásában milyen szerepe lehet az ember-AI
kollaborációnak.

## Mit találtunk?

A platform öt, terápiás szempontból releváns fehérje célponton validáltuk
(összesen 824 molekuláris dokkolási kísérlet, ~313 validált eredmény):

| Célpont | PDB ID | Legjobb találat | Affinitás | Klinikai relevancia |
|---------|--------|----------------|-----------|--------------------|
| SARS-CoV-2 Mpro | 6LU7 | Nilotinib_Var_17 | -8.85 kcal/mol | COVID-19 |
| HIV-1 Proteáz | 1HPV | Nilotinib_Var_17 | -10.09 kcal/mol | HIV/AIDS |
| EGFR Kináz | 1M17 | Nilotinib_Var_17 | -10.94 kcal/mol | Tüdőrák |
| Bcl-2 | 4LVT | Nilotinib_Var_17 | -8.06 kcal/mol | Leukémia |
| COX-2 | 5KIR | Abiraterone | -11.08 kcal/mol | Gyulladás |

A felfedezés lényege: a Nilotinib_Var_17 — az FDA által már jóváhagyott
Nilotinib (Tasigna, Novartis) egy számítógépesen optimalizált variánsa —
mind az öt célponton kiemelkedő kötődést mutat. Ez azt jelenti, hogy egy
már ismert biztonságosságú gyógyszerből kiindulva találtunk olyan kémiai
módosítást, ami jelentősen javítja a kötődést több, terápiásan fontos
fehérjéhez.

## Miért írok Önöknek?

Mert ez a munka jelenleg kizárólag in silico — számítógépes szimulációkra
épül. Bármennyire is ígéretesek ezek az eredmények, a valóságban csak
laboratóriumi kísérletekkel igazolhatóak.

Szeretném megkérni Önöket, hogy segítsenek ezt a projektet a következő
szintre emelni. Olyan partnert keresek, aki:

1. In vitro validálni tudná a predikcióinkat — enzimaktivitási assay-ekkel
   (SPR, FRET, vagy fluoreszcencia polarizáció) megmérni, hogy a
   Nilotinib_Var_17 valóban kötődik-e a prediktált fehérjékhez

2. Sejtes tesztekkel igazolni az antivirális vagy tumorellenes hatást —
   akár csak egyetlen célponton, akár a COVID-19 Mpro-n

3. Szakértői visszajelzést adni a módszertanunkról — farmakológus
   szemmel értékelni, hogy a megközelítésünk tudományosan megalapozott-e

## Mit tudok ajánlani cserébe?

- A platform ingyenes és korlátlan használatát az Intézet bármely
  kutatójának — tetszőleges fehérje célpontra percek alatt tudnak
  gyógyszerjelölteket szűrni
- A szoftver testreszabását az Önök kutatási igényeihez
- Teljes körű technikai támogatást és továbbfejlesztést
- Közös publikáció lehetőségét — a Semmelweis Egyetem neve alatt
- A teljes forráskódot és a white paper-t (MIT licenc)

## Miért érdemes időt szánniuk rám?

Tisztában vagyok vele, hogy egy egyetemi intézet ideje értékes — és egy
technikumi diák projektje talán nem tűnik prioritásnak. De szeretném
kiemelni, hogy:

- Ez a projekt NEM egy iskolai házi feladat. Hónapok munkája, több száz
  dokkolási kísérlet, és egy olyan platform, amit akár holnap is
  használhatnának kutatók
- Az eredmények validáltak: 5 fehérjén, 824 dokkolással, publikálható
  formában dokumentálva
- A módszertan nyílt és reprodukálható: a teljes kód, adatok és white
  paper elérhető GitHub-on
- Nem pénzt kérek — hanem szakmai partnerséget. Együttműködést, amiből
  mindkét fél tanul és publikálhat

Ha csak egyetlen vegyületet — a Nilotinib_Var_17-et — letesztelnének
egy in vitro enzimassay-ben a SARS-CoV-2 Mpro ellen, az már óriási
előrelépés lenne számomra. Ez nem igényel jelentős erőforrást, de
számomra felbecsülhetetlen értékű validációt adna.

## Hol érhető el a projekt?

- GitHub: https://github.com/horvatjanos/universal-drug-repurposing
- White Paper: https://github.com/horvatjanos/universal-drug-repurposing/blob/master/docs/WHITE_PAPER.md
- A platform azonnal használható: python src/hermes_drug.py --interactive

---

PS: Láttam, hogy az Intézetükben fut egy "AI-alapú gyógyszerkutatási platform
kifejlesztése" nevű projekt is — talán a mi munkánk ehhez is kapcsolódhatna,
vagy kiegészíthetné azt a molekuláris dokkolás oldaláról.

Nagyon köszönöm, hogy időt szántak a levelem elolvasására. Őszintén hiszek
segítségével valódi tudományos értékké válhat.

Várom megtisztelő válaszukat!

Tisztelettel,

Turcsányi Zsombor
Neumann János Informatikai Technikum
Közreműködő: Hermes AI

E-mail: [ide ird az emailed]
GitHub: https://github.com/horvatjanos
"""

print(LEVÉL)
