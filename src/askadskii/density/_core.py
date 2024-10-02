import math
import sys
from collections import defaultdict
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

from askadskii.logging import DEBUG, get_child_logger
from askadskii.utils.dataset import decrypt_manager

SMILES: TypeAlias = str

K_CONSTANTS = 0.681
UNIT_CONSTANTS = (
    scipy.constants.Avogadro
    / (scipy.constants.centi / scipy.constants.angstrom) ** 3
)

_logger = get_child_logger(__name__)


COL_VDW_VOLUME = "V_i"

DIRPATH_DATA = Path(__file__).parent
FILEPATH_DISTANCE = DIRPATH_DATA / "distance.csv"
FILEPATH_RADII = DIRPATH_DATA / "radii.csv"
FILEPATH_COVALENT = DIRPATH_DATA / "radius-covalent.csv"
# FILEPATH_DATA = DIRPATH_DATA / "askadskii.csv"
# if not FILEPATH_DATA.is_file():
#     FILEPATH_DATA = DIRPATH_DATA / "askadskii_hydrocarbon.csv"


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
    use_estimated_params: bool = False,
) -> float:
    """estimate van der waals' volume of polymers

    Parameters
    ----------
    smiles_or_mol : Union[SMILES, rdkit.Chem.rdchem.Mol]
        SMILES or rdkit.Chem.rdchem.Mol
    unit : Literal["cm^3/mol", "angstrom^3"], optional
        unit, by default "cm^3/mol"
    use_estimated_params : bool, optional
        use estimated covalent bond dinstance or not, by default False

    Returns
    -------
    float
        estimated van der waals' volume

    Examples
    --------
    >>> from askadskii.density import estimate_vdw_volume
    >>> estimate_vdw_volume("*CC*")

    """
    mol = _get_mol(smiles_or_mol)
    map_head_tail = _get_head_tail(mol)
    mol_with_hs = Chem.AddHs(mol)

    # idxがAddHsによって変わっていないか念の為確認
    assert all(
        _atom.GetAtomMapNum() - 1 == _atom.GetIdx()
        for _atom in mol_with_hs.GetAtoms()
        if _atom.GetSymbol() != "H"
    )

    with decrypt_manager():
        # 結合距離情報
        df_distance = pd.read_csv(
            FILEPATH_DISTANCE, encoding="utf-8", comment="#"
        )
        df_distance.dropna(axis=0, how="all", inplace=True)
        map_distance = MappingProxyType(
            {
                _get_key_distance(
                    _sr_row["Atom_1"],
                    _sr_row["Atom_2"],
                    _sr_row["BondType"],
                ): _sr_row["Distance"]
                for _, _sr_row in df_distance.iterrows()
            }
        )

        # 半径情報
        df_radii = pd.read_csv(
            FILEPATH_RADII, index_col=0, encoding="utf-8", comment="#"
        )
        map_radii = MappingProxyType(df_radii.iloc[:, 0].to_dict())

        # 結合距離が定義されていない場合
        if use_estimated_params:
            df_radii_coval = pd.read_csv(
                FILEPATH_COVALENT, index_col=0, encoding="utf-8", comment="#"
            )
            df_radii_coval.loc[:, "1.5"] = (
                df_radii_coval.loc[:, ["1.0", "2.0"]]
                .mean(axis=1)
                .rename("1.5")
            )
            map_radii_coval = MappingProxyType(
                df_radii_coval.to_dict(orient="index")
            )

    vdw_volume = 0.0
    for atom in mol_with_hs.GetAtoms():
        atom: Chem.rdchem.Atom
        symbol = atom.GetSymbol()

        # 繰り返し末端は無視
        if symbol == "*":
            continue

        if atom.GetIsAromatic():
            symbol = symbol.lower()

        radii = map_radii[symbol.capitalize()]
        _vdw_volume_i = (4 / 3) * math.pi * (radii**3)

        for bond in atom.GetBonds():
            bond: Chem.rdchem.Bond
            if bond.GetBeginAtomIdx() == atom.GetIdx():
                atom_other = bond.GetEndAtom()
            elif bond.GetEndAtomIdx() == atom.GetIdx():
                atom_other = bond.GetBeginAtom()
            else:
                raise ValueError

            symbol_other = atom_other.GetSymbol()
            # 繰り返し末端に隣接している場合
            if symbol_other == "*":
                # もう1つの繰り返し末端の隣の原子に置き換える
                atom_other: Chem.rdchem.Atom = mol_with_hs.GetAtomWithIdx(
                    map_head_tail[atom_other.GetIdx()]
                ).GetNeighbors()[0]  # 1つしか結合していないのは確認済み
                symbol_other = atom_other.GetSymbol()

            if symbol != "H" and atom_other.GetIsAromatic():
                symbol_other = symbol_other.lower()
            radii_other = map_radii[symbol_other.capitalize()]

            _key_bond = _get_key_distance(
                symbol, symbol_other, bond.GetBondTypeAsDouble()
            )
            if _key_bond in set(map_distance.keys()):
                distance = map_distance[_key_bond]
            elif use_estimated_params:
                _logger.warning(f"{_key_bond} is estiamted.")
                distance = sum(
                    map_radii_coval[_symbol.capitalize()][
                        f"{bond.GetBondTypeAsDouble():.1f}"
                    ]
                    for _symbol in (symbol, symbol_other)
                )
            else:
                raise KeyError(
                    f"{_key_bond} is not defined; check 'use_estimated_params' mode."
                )

            h_i = radii - (radii**2 + distance**2 - radii_other**2) / (
                2 * distance
            )
            _vdw_volume_i -= (1 / 3) * math.pi * (h_i**2) * (3 * radii - h_i)

        # nm to angstrom
        _vdw_volume_i *= (scipy.constants.nano / scipy.constants.angstrom) ** 3
        # _vdw_volume_i = float(f"{_vdw_volume_i:.1f}")
        _logger.debug(
            f"{symbol} {tuple(_atom.GetSymbol() for _atom in atom.GetNeighbors())}, {_vdw_volume_i:.1f}"
        )

        vdw_volume += _vdw_volume_i

    if unit == "cm^3/mol":
        vdw_volume *= UNIT_CONSTANTS
    elif unit == "angstrom^3":
        pass
    else:
        raise ValueError("'cm^3/mol' or 'angstrom^3'")
    return vdw_volume


def estimate_density(
    smiles_or_mol: Union[SMILES, Chem.rdchem.Mol],
    use_estimated_params: bool = False,
) -> float:
    """estimate density of polymers

    Parameters
    ----------
    smiles_or_mol : Union[SMILES, Chem.rdchem.Mol]
        SMILES or rdkit.Chem.rdchem.Mol
    use_estimated_params : bool, optional
        use estimated covalent bond dinstance or not, by default False

    Returns
    -------
    float
        estimated polymer's density

    Examples
    --------
    >>> from askadskii.density import estimate_density
    >>> estimate_density("*CC*")
    """
    mol = _get_mol(smiles_or_mol)
    vdw_volume = estimate_vdw_volume(
        mol, unit="cm^3/mol", use_estimated_params=use_estimated_params
    )

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


def _get_key_distance(
    symbol1: str, symbol2: str, bondtype: float
) -> tuple[str, str, str]:
    if symbol1 == "H":
        symbol2 = symbol2.capitalize()
    elif symbol2 == "H":
        symbol1 = symbol1.capitalize()

    return tuple(sorted((symbol1, symbol2)) + ["{0:.1f}".format(bondtype)])


if __name__ == "__main__":
    _logger.setLevel(DEBUG)

    # for _smiles in ("*CC(*)(C)C", "*CC(C#N)*", "*CC(c1ccccc1)*"):
    for _smiles in (
        # "*CC(*)(C)C",
        # "*CC(C#N)*",
        # "*CC(*)(C)C(=O)OC",
        # "*CC(*)(C)C(=O)OCC",
        "*CC(c1ccccc1)*",
        # "*CC(c1ccccc1)(C)*",
        # "C(*)C(*)(C)C(=O)OC",
        # "CC(*)C(*)(C)C(=O)OC",
        # "*CC*",
        # "*NCCCCCCNC(=O)CCCCC(=O)*",
        # "[1*]SCC(CS[2*])(CSC(=S)Nc1ccc(NC([2*])=S)cc1)CSC(=S)Nc2ccc(NC([1*])=S)cc2",
        # "*CC(O)*",
    ):
        _logger.debug(_smiles)
        # print("{:.1f}".format(estimate_vdw_volume(_smiles, unit="cm^3/mol")))
        # print("{:.1f}".format(estimate_vdw_volume(_smiles, unit="angstrom^3")))
        print("{:.2f}".format(estimate_density(_smiles)))
