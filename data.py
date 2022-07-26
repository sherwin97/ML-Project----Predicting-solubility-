import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

data = pd.read_csv(
    "/Users/sherwinng/Desktop/VS_code/chem/ci034243xsi20040112_053635.txt"
)

def desc(smiles):
    """
    Take in a list of rdkit object
    Allow user to obtain LogP, molecular weight and no of rotational
    bonds
    """
    mol_list = [Chem.MolFromSmiles(mol) for mol in data.SMILES]

    base_data = np.arange(1, 1)
    i = 0
    for mol in mol_list:
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_Molwt = Descriptors.MolWt(mol)
        desc_NumRotableBonds = Descriptors.NumRotatableBonds(mol)
        desc_TPSA = Descriptors.TPSA(mol)

        row = np.array([desc_MolLogP, desc_Molwt, desc_NumRotableBonds, desc_TPSA])

        if i == 0:
            base_data = row

        else:
            base_data = np.vstack([base_data, row])

        i += 1

    column_names = ["MolLogP", "MolWt", "NumRotableBonds", "TPSA"]
    descriptors = pd.DataFrame(data=base_data, columns=column_names)

    return descriptors


def AromaticAtoms(m):
    """
    Take in an rdkit object, obtain the number of aromatic atoms in a molecule
    """
    aromatic_atoms2 = [
        m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())
    ]
    aa_count = []
    for i in aromatic_atoms2:
        if i == True:
            aa_count.append(i)

    sum_aa_count = sum(aa_count)
    return sum_aa_count


mol_list = [Chem.MolFromSmiles(mol) for mol in data.SMILES]
df = desc(data.SMILES)  # obtain LogP, Molwt and no. of rot bonds
desc_AromaticProportion = [
    AromaticAtoms(mol) / Descriptors.HeavyAtomCount(mol) for mol in mol_list
]  # obtain aromatic proportion
df_desc_AromaticProportion = pd.DataFrame(
    desc_AromaticProportion, columns=["Aromatic Proportion"]
)


X = pd.concat([df, df_desc_AromaticProportion], axis=1)
y = data.iloc[:, 1]

# to check on the correlation between each factor and solubility. For this project, as a rule of thumb, more than .5 = correlated.
overall_data = pd.concat([df, df_desc_AromaticProportion, data.iloc[:, 1]], axis=1)
