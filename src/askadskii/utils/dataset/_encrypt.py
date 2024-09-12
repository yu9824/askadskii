import glob
import importlib
from logging import DEBUG
from pathlib import Path
from typing import Optional, Union

from cryptography.fernet import Fernet

from askadskii.logging import get_child_logger
from askadskii.logging._logging import _get_root_logger_name
from askadskii.utils.dataset._constants import KEY_SPLIT

DIRPATH_ROOT = Path(
    importlib.import_module(_get_root_logger_name()).__file__
).parent

_logger = get_child_logger(__name__)


def __encrypt(key: Optional[Union[bytes, str]] = None) -> bytes:
    if key is None:
        key = Fernet.generate_key()

    contents_all = bytes()
    for filepath_csv in glob.glob(str(DIRPATH_ROOT / "**" / "*.csv")):
        filepath_csv = Path(filepath_csv)

        with open(filepath_csv, mode="rb") as f:
            _content = f.read()
        contents_all += (
            bytes(filepath_csv.relative_to(DIRPATH_ROOT))
            + KEY_SPLIT
            + _content
            + KEY_SPLIT
        )

    cipher_suite = Fernet(key)
    contents_all_crypted = cipher_suite.encrypt(contents_all)

    with open(Path(__file__).parent / ".dataset", mode="wb") as f:
        f.write(contents_all_crypted)

    _logger.debug("key: '{}'".format(key.decode("utf-8")))
    return key


if __name__ == "__main__":
    _logger.setLevel(DEBUG)
    __encrypt()
