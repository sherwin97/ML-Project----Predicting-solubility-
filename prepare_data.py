import argparse

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split


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


def prepare_data(
    smiles,
    solubility,
    output_path1,
    output_path2,
    output_path3,
    output_path4,
    output_path5,
):
    """
    Load csv containing smiles. Convert str to rdkit object and obtain MolLogP, MolWt, NumRotBonds, TPSA and AromaticProp. Return a dataframe.
    """
    smiles_list = [item for item in open(smiles).read().replace("\n", ",").split(",")]
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

    y = pd.read_csv(solubility)

    X_train, X_test, y_train, y_test = train_test_split(
        df, y, test_size=0.2, random_state=123
    )

    return (
        df.to_csv(output_path1),
        X_train.to_csv(output_path2, index = False),
        X_test.to_csv(output_path3, index = False),
        y_train.to_csv(output_path4, index = False),
        y_test.to_csv(output_path5, index=False),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--X_input", help="Enter the file path containing SMILES")
    parser.add_argument("--y_input", help="Enter the file path containing solubilities")
    parser.add_argument(
        "--output",
        help="Enter the file path to save csv containing molecular descriptions",
    )
    parser.add_argument("--X_train", help="Enter the file path to save X_train")
    parser.add_argument("--X_test", help="Enter the file path to save X_test")
    parser.add_argument("--y_train", help="Enter the file path to save y_train")
    parser.add_argument("--y_test", help="Enter the file path to save y_test")
    args = parser.parse_args()

    prepare_data(
        args.X_input,
        args.y_input,
        args.output,
        args.X_train,
        args.X_test,
        args.y_train,
        args.y_test,
    )
