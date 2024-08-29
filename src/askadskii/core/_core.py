from collections import Counter, defaultdict
from pathlib import Path
from types import MappingProxyType
from typing import TypeAlias, Union

import pandas as pd
import rdkit.Chem
import scipy.constants
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

SMILES: TypeAlias = str

K_CONSTANTS = 0.681
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
    ).dropna(axis=0, how="any")

    df_condition = df_data.drop(COL_VDW_VOLUME, axis=1).astype(
        {
            _col: int
            for _col in df_data.columns
            if _col not in {COL_VDW_VOLUME, "Symbol"}
        }
    )
    sr_v = df_data.loc[:, COL_VDW_VOLUME]

    vdw_volume = 0.0
    for atom in mol_with_hs.GetAtoms():
        atom: Chem.rdchem.Atom
        symbol = atom.GetSymbol()

        # 繰り返し末端は無視
        if symbol == "*":
            continue

        cnt = Counter()
        for atom_neighbor in atom.GetNeighbors():
            atom_neighbor: Chem.rdchem.Atom

            symbol_neighbor = atom_neighbor.GetSymbol()
            # 繰り返し末端に隣接している場合
            if symbol_neighbor == "*":
                # もう1つの繰り返し末端の隣の原子に置き換える
                atom_neighbor = mol_with_hs.GetAtomWithIdx(
                    map_head_tail[atom_neighbor.GetIdx()]
                ).GetNeighbors()[0]  # 1つしか結合していないのは確認済み
                symbol_neighbor = atom_neighbor.GetSymbol()

            if symbol != "H" and atom_neighbor.GetIsAromatic():
                symbol_neighbor = symbol_neighbor.lower()
            cnt[symbol_neighbor] += 1

        dict_info = dict(cnt)
        dict_info["Symbol"] = symbol

        for _key in set(df_condition.columns) - set(dict_info.keys()):
            dict_info[_key] = 0
        sr_info = pd.Series(dict_info).loc[df_condition.columns]

        mask = df_condition.eq(sr_info, axis=1).all(axis=1)
        assert mask.sum() == 1

        vdw_volume += sr_v[mask].item()
    return vdw_volume


def estimate_density(smiles_or_mol) -> float:
    mol = _get_mol(smiles_or_mol)
    vdw_volume = estimate_vdw_volume(mol)

    return (K_CONSTANTS * MolWt(mol)) / (
        scipy.constants.Avogadro
        * vdw_volume
        * (scipy.constants.angstrom / scipy.constants.centi) ** 3
    )


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
    print(estimate_density("[1*]CC(C[2*])(C[2*])C[1*]"))
    print(estimate_density("[*]CC[*]"))
    print(estimate_density("CCC"))

    # _get_head_tail("[*]CC[*]")
