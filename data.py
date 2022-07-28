import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

# load data set, data obtained from delaney paper on ESOL
data = pd.read_csv(
    "ci034243xsi20040112_053635.txt"
)

# to obtain LogP, MW, RB TPSA for each mol


def desc(smiles):
    """
    Take in a list of rdkit object
    Allow user to obtain LogP, molecular weight, no of rotational bonds, TPSA and return a df
    """
    mol_list = [Chem.MolFromSmiles(mol) for mol in data.SMILES]

    mol_MolLogP_list = [Descriptors.MolLogP(mol) for mol in mol_list]
    mol_MolWt_list = [Descriptors.MolWt(mol) for mol in mol_list]
    mol_NumRotableBonds_list = [Descriptors.NumRotatableBonds(mol) for mol in mol_list]
    mol_TPSA_list = [Descriptors.TPSA(mol) for mol in mol_list]

    df = pd.DataFrame(
        {
            "MolLogP": mol_MolLogP_list,
            "MolWt": mol_MolWt_list,
            "NumRotableBonds": mol_NumRotableBonds_list,
            "TPSA": mol_TPSA_list,
        }
    )

    return df


def AromaticAtoms(m):
    """
    Take in an rdkit object, obtain the number of aromatic atoms in a molecule
    Aromatic atoms per molecule are obtain and returned. Only heavy atoms are considered.
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
