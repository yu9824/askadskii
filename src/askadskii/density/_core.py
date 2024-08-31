import sys
from collections import Counter, defaultdict
from pathlib import Path
from types import MappingProxyType
from typing import TypeAlias, Union

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

import pandas as pd
import rdkit.Chem
import scipy.constants
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

SMILES: TypeAlias = str

K_CONSTANTS = 0.681
UNIT_CONSTANTS = (
    scipy.constants.Avogadro
    / (scipy.constants.centi / scipy.constants.angstrom) ** 3
)


COL_VDW_VOLUME = "V_i"

DIRPATH_DATA = Path(__file__).parent
FILEPATH_DATA = DIRPATH_DATA / "askadskii.csv"
if not FILEPATH_DATA.is_file():
    FILEPATH_DATA = DIRPATH_DATA / "askadskii_hydrocarbon.csv"


def _get_mol(
    smiles_or_mol: Union[SMILES, rdkit.Chem.rdchem.Mol],
) -> rdkit.Chem.rdchem.Mol:
    if isinstance(smiles_or_mol, rdkit.Chem.rdchem.Mol):
        mol = smiles_or_mol
    elif isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        raise ValueError

    for map_num, atom in enumerate(
        mol.GetAtoms(),
        start=1,  # デフォルトは0のため、重複を避ける
    ):
        atom: Chem.rdchem.Atom
        atom.SetAtomMapNum(map_num)

    return mol


def estimate_vdw_volume(
    smiles_or_mol: Union[SMILES, rdkit.Chem.rdchem.Mol],
    unit: Literal["cm^3/mol", "angstrom^3"] = "cm^3/mol",
) -> float:
    mol = _get_mol(smiles_or_mol)
    map_head_tail = _get_head_tail(mol)
    mol_with_hs = Chem.AddHs(mol)

    # idxがAddHsによって変わっていないか念の為確認
    assert all(
        _atom.GetAtomMapNum() - 1 == _atom.GetIdx()
        for _atom in mol_with_hs.GetAtoms()
        if _atom.GetSymbol() != "H"
    )

    df_data = pd.read_csv(
        FILEPATH_DATA, index_col=0, delimiter=",", encoding="utf-8"
    )
    df_data = df_data.loc[~df_data.loc[:, COL_VDW_VOLUME].isnull()]

    df_condition = df_data.drop([COL_VDW_VOLUME, "comment"], axis=1).fillna(
        0.0
    )
    sr_v = df_data.loc[:, COL_VDW_VOLUME]

    vdw_volume = 0.0
    for atom in mol_with_hs.GetAtoms():
        atom: Chem.rdchem.Atom
        symbol = atom.GetSymbol()

        # 繰り返し末端は無視
        if symbol == "*":
            continue

        if atom.GetIsAromatic():
            symbol = symbol.lower()

        cnt = Counter()
        for bond in atom.GetBonds():
            bond: Chem.rdchem.Bond
            if bond.GetBeginAtomIdx() == atom.GetIdx():
                atom_another = bond.GetEndAtom()
            elif bond.GetEndAtomIdx() == atom.GetIdx():
                atom_another = bond.GetBeginAtom()
            else:
                raise ValueError

            symbol_another = atom_another.GetSymbol()
            # 繰り返し末端に隣接している場合
            if symbol_another == "*":
                # もう1つの繰り返し末端の隣の原子に置き換える
                atom_another = mol_with_hs.GetAtomWithIdx(
                    map_head_tail[atom_another.GetIdx()]
                ).GetNeighbors()[0]  # 1つしか結合していないのは確認済み
                symbol_another = atom_another.GetSymbol()

            if symbol != "H" and atom_another.GetIsAromatic():
                symbol_another = symbol_another.lower()
            cnt[symbol_another] += bond.GetBondTypeAsDouble()
        dict_info = dict(cnt)
        dict_info["Symbol"] = symbol
        dict_info["N_Bonds"] = len(atom.GetBonds())

        for _key in set(df_condition.columns) - set(dict_info.keys()):
            dict_info[_key] = 0
        sr_info = pd.Series(dict_info).loc[df_condition.columns]

        mask = df_condition.eq(sr_info, axis=1).all(axis=1)
        assert mask.sum() == 1, f"{sr_info=}"

        vdw_volume += sr_v[mask].item()
    if unit == "cm^3/mol":
        vdw_volume *= UNIT_CONSTANTS
    elif unit == "angstrom^3":
        pass
    else:
        raise ValueError("'cm^3/mol' or 'angstrom^3'")
    return vdw_volume


def estimate_density(smiles_or_mol) -> float:
    mol = _get_mol(smiles_or_mol)
    vdw_volume = estimate_vdw_volume(mol, unit="cm^3/mol")

    return (K_CONSTANTS * MolWt(mol)) / vdw_volume


def _get_head_tail(
    smiles_or_mol: Union[SMILES, rdkit.Chem.rdchem.Mol],
) -> "MappingProxyType[int, int]":
    _dict_pairs = defaultdict(list)
    mol = _get_mol(smiles_or_mol)
    for atom in mol.GetAtoms():
        atom: Chem.rdchem.Atom
        if atom.GetSymbol() == "*":
            # 隣接原子が1つしかないことを確認
            assert len(atom.GetNeighbors()) == 1

            _dict_pairs[atom.GetIsotope()].append(
                atom.GetAtomMapNum() - 1
            )  # indexは0-index
    # 2つで1ペアになっていることを確認する
    assert all(len(_pair) == 2 for _pair in _dict_pairs.values())

    dict_head_tail = dict()
    for head, tail in _dict_pairs.values():
        dict_head_tail[head] = tail
        dict_head_tail[tail] = head
    return MappingProxyType(dict_head_tail)


if __name__ == "__main__":
    print(estimate_vdw_volume("*CC(C)=CC*"))
    print(estimate_density("*CC(C)=CC*"))

    print(estimate_vdw_volume("*CC(C#N)*"))
    print(estimate_density("*CC(C#N)*"))

    print(estimate_vdw_volume("*CC(c1ccccc1)*"))
    print(estimate_density("*CC(c1ccccc1)*"))
