import random
from rdkit import Chem
from rdkit.Chem import Descriptors

class LabRat:
    def __init__(self):
        # Simple dictionary to map common drug names to SMILES for demo purposes
        self.molecule_db = {
            "ciprofloxacin": "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
            "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
            "chimera-x1": "C1=NC2=C(N1)C(=NC=N2)N", # Mock SMILES
            "penicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
        }

    def run_simulation(self, molecule_name):
        """
        Simulates a drug screening process.
        1. Resolves molecule name to SMILES.
        2. Calculates properties using RDKit.
        3. Generates a mock binding affinity score.
        """
        molecule_name = molecule_name.lower().strip()
        smiles = self.molecule_db.get(molecule_name)

        if not smiles:
            # Try to treat input as SMILES if not found in DB
            mol = Chem.MolFromSmiles(molecule_name)
            if mol:
                smiles = molecule_name
            else:
                return f"Lab Rat Error: Molecule '{molecule_name}' not found in database and invalid as SMILES."

        # Calculate Properties using RDKit
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        # Mock Binding Affinity (Randomized for demo, but deterministic for specific molecules)
        # In a real app, this would call AutoDock Vina
        random.seed(molecule_name) 
        binding_score = round(random.uniform(0.5, 0.99), 2)
        
        # Determine "Success" based on score
        status = "High Affinity" if binding_score > 0.8 else "Moderate Affinity" if binding_score > 0.6 else "Low Affinity"

        report = f"""
🧪 **Lab Rat Simulation Report**
**Molecule:** {molecule_name.capitalize()}
**SMILES:** `{smiles}`

**Molecular Properties (RDKit):**
- Molecular Weight: {mw:.2f} g/mol
- LogP (Lipophilicity): {logp:.2f}
- H-Bond Donors: {hbd}
- H-Bond Acceptors: {hba}

**Screening Results:**
- **Binding Affinity Score:** {binding_score} ({status})
- **Lipinski's Rule of 5:** {"Pass" if logp <= 5 and mw <= 500 and hbd <= 5 and hba <= 10 else "Fail"}
"""
        return report

# Test
if __name__ == "__main__":
    rat = LabRat()
    print(rat.run_simulation("ciprofloxacin"))
