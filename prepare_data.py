import argparse

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors



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
        if i:
            aa_count.append(i)

    sum_aa_count = sum(aa_count)
    return sum_aa_count


def prepare_data(
    path_smiles,
    path_features
):
    """
    Load csv containing smiles. Convert str to rdkit object and obtain MolLogP, MolWt, NumRotBonds, TPSA and AromaticProp. Return a dataframe.
    """
    smiles_list = [item for item in open(path_smiles).read().replace("\n", ",").split(",")]
    mol_list = [Chem.MolFromSmiles(mol) for mol in smiles_list]

    mol_MolLogP_list = [Descriptors.MolLogP(mol) for mol in mol_list]
    mol_MolWt_list = [Descriptors.MolWt(mol) for mol in mol_list]
    mol_NumRotableBonds_list = [Descriptors.NumRotatableBonds(mol) for mol in mol_list]
    mol_TPSA_list = [Descriptors.TPSA(mol) for mol in mol_list]
    mol_AromaticProportion = [
        AromaticAtoms(mol) / Descriptors.HeavyAtomCount(mol) for mol in mol_list
    ]
    df = pd.DataFrame(
        {
            "MolLogP": mol_MolLogP_list,
            "MolWt": mol_MolWt_list,
            "NumRotableBonds": mol_NumRotableBonds_list,
            "TPSA": mol_TPSA_list,
            "AromaticProportion": mol_AromaticProportion,
        }
    )

    return (
        df.to_csv(path_features, header=False, index=False)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_smiles", help="Enter the file path containing SMILES")
    parser.add_argument(
        "--path_features",
        help="Enter the file path to save csv containing molecular descriptions")

    args = parser.parse_args()

    prepare_data(args.path_smiles, args.path_features)
