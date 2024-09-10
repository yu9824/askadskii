import glob
import importlib
from pathlib import Path
from typing import Union

from cryptography.fernet import Fernet

from askadskii.logging._logging import _get_root_logger_name
from askadskii.utils.dataset._constants import KEY_SPLIT

DIRPATH_ROOT = Path(
    importlib.import_module(_get_root_logger_name()).__file__
).parent


def __decrypt(key: Union[str, bytes]) -> None:
    with open(Path(__file__).parent / "dataset", mode="rb") as f:
        contents_raw = f.read()

    cipher_suite = Fernet(key)
    contents = cipher_suite.decrypt(contents_raw)

    contents_split = contents.split(KEY_SPLIT)[:-1]
    for i_file in range(contents.count(KEY_SPLIT) // 2):
        with open(contents_split[2 * i_file], mode="wb") as f:
            f.write(contents_split[2 * i_file + 1])


if len(glob.glob(str(DIRPATH_ROOT / "**" / "*.csv"))) == 0:
    __decrypt(input("key: "))
