import contextlib
import glob
import importlib
import os
from pathlib import Path
from typing import Optional, Union

from cryptography.fernet import Fernet
from dotenv import load_dotenv

from askadskii.logging._logging import _get_root_logger_name
from askadskii.utils.dataset._constants import KEY_ASKADSKII, KEY_SPLIT

DIRPATH_ROOT = Path(
    importlib.import_module(_get_root_logger_name()).__file__
).parent
FILEPATH_ENV = DIRPATH_ROOT / ".env"


def __decrypt(key: Union[str, bytes]) -> None:
    with open(Path(__file__).parent / ".dataset", mode="rb") as f:
        contents_raw = f.read()

    cipher_suite = Fernet(key)
    contents = cipher_suite.decrypt(contents_raw)

    contents_split = contents.split(KEY_SPLIT)[:-1]
    for i_file in range(contents.count(KEY_SPLIT) // 2):
        with open(
            DIRPATH_ROOT / contents_split[2 * i_file].decode("utf-8"),
            mode="wb",
        ) as f:
            f.write(contents_split[2 * i_file + 1])


@contextlib.contextmanager
def decrypt_manager(key: Optional[Union[str, bytes]] = None):
    try:
        if FILEPATH_ENV.is_file():
            load_dotenv(FILEPATH_ENV)

        if key:
            __decrypt(key)
        elif KEY_ASKADSKII in set(os.environ):
            __decrypt(os.environ[KEY_ASKADSKII])
        else:
            _activation_key = input("activation key: ")
            __decrypt(_activation_key)

            with open(FILEPATH_ENV, mode="w", encoding="utf-8") as f:
                f.write(f"{KEY_ASKADSKII}={_activation_key}")

        yield

    finally:
        for _filepath_dataset in glob.glob(str(DIRPATH_ROOT / "**" / "*.csv")):
            os.remove(_filepath_dataset)


if __name__ == "__main__":
    if FILEPATH_ENV.is_file():
        load_dotenv(FILEPATH_ENV)

    if KEY_ASKADSKII in set(os.environ):
        __decrypt(os.environ[KEY_ASKADSKII])
    else:
        _activation_key = input("activation key: ")
        __decrypt(_activation_key)

        with open(FILEPATH_ENV, mode="w", encoding="utf-8") as f:
            f.write(f"{KEY_ASKADSKII}={_activation_key}")
