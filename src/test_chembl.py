from chembl_webresource_client.new_client import new_client
import sys

def test_chembl_connection():
    print("=== ChEMBL API Kapcsolat Teszt ===")
    print("Csatlakozás folyamatban...")
    try:
        molecule = new_client.molecule
        # Megpróbálunk lekérni 5 darab FDA-engedélyezett gyógyszert tesztnek
        print("Adatok lekérése a szerverről...")
        res = molecule.filter(max_phase=4).only(['molecule_chembl_id', 'pref_name'])[:5]
        
        found = False
        for i, mol in enumerate(res):
            found = True
            name = mol.get('pref_name') or "Névtelen"
            print(f"{i+1}. Talált gyógyszer: {name} (ID: {mol['molecule_chembl_id']})")
            
        if found:
            print("\nSIKER! A szerver válaszol, a kapcsolat stabil.")
        else:
            print("\nFIGYELEM: A szerver válaszolt, de nem küldött adatot.")
            
        return True
    except Exception as e:
        print(f"\nHIBA történt a csatlakozás során!")
        print(f"Részletek: {str(e)}")
        if "500" in str(e):
            print("\nTipp: Ez egy szerveroldali hiba (Internal Server Error). A ChEMBL szerverei valószínűleg túlterheltek. Érdemes 10-20 percet várni.")
        return False

if __name__ == "__main__":
    test_chembl_connection()
