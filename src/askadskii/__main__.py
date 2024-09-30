import argparse
import importlib
import sys
from logging import DEBUG, getLogger
from pathlib import Path
from typing import Optional

if sys.version_info >= (3, 9):
    from collections.abc import Sequence
else:
    from typing import Sequence


from askadskii import __version__
from askadskii.density import estimate_density, estimate_vdw_volume
from askadskii.logging._logging import _get_root_logger_name

__all__ = ("main",)

root_logger = getLogger(_get_root_logger_name())


def main(cli_args: Sequence[str], prog: Optional[str] = None) -> None:
    parser = argparse.ArgumentParser(prog=prog, description="")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        help="show current version",
        version=f"%(prog)s: {__version__}",
    )

    subparsers = parser.add_subparsers()
    parser_density = subparsers.add_parser(
        "density", help="see 'density --help'"
    )
    parser_vdw_volume = subparsers.add_parser(
        "vdw-volume",
        help="see 'vdw-volume --help'",
        description="Unit is 'cm^3/mol'",
    )
    parser_density.set_defaults(func=estimate_density)
    parser_vdw_volume.set_defaults(func=estimate_vdw_volume)
    for _subparser in (parser_density, parser_vdw_volume):
        _subparser.add_argument(
            "smiles", type=str, nargs="*", help="smiles", metavar="SMILES"
        )
        _subparser.add_argument(
            "-i",
            "--input-file",
            type=Path,
            help="A filepath containing a list of SMILES separated by newlines",
        )
        _subparser.add_argument(
            "-d",
            "--decimal",
            type=int,
            default=3,
            help="the number of decimal places to display, by default 3",
        )
        _subparser.add_argument(
            "--debug", action="store_true", help="debug mode"
        )
        _subparser.add_argument(
            "--deprecated", action="store_true", help="use deprecated method"
        )
        _subparser.add_argument(
            "-o", "--output-file", type=Path, help="filepath output"
        )
        _subparser.add_argument(
            "-e",
            "--use-estimated-params",
            action="store_true",
            help="use estimated covalent bond dinstance",
        )

    _common(parser.parse_args(cli_args))


def entrypoint() -> None:
    main(sys.argv[1:])


def _common(args: argparse.Namespace) -> None:
    list_smiles: "list[str]" = args.smiles
    filepath_input: Optional[Path] = args.input_file
    filepath_output: Optional[Path] = args.output_file

    if args.deprecated:
        _deprecated_module = importlib.import_module(
            f"{_get_root_logger_name()}.density._deprecated"
        )
        func = getattr(_deprecated_module, args.func.__name__)
    else:
        func = args.func

    if args.debug:
        root_logger.setLevel(DEBUG)

    if not (bool(list_smiles) or bool(filepath_input)):
        raise ValueError("one of them should be specified.")

    elif bool(list_smiles) and bool(filepath_input):
        raise ValueError("Either of them should be specified.")

    elif filepath_input:
        with open(filepath_input, mode="r", encoding="utf-8") as f:
            tup_smiles = tuple(f.read().split())

    else:
        tup_smiles = tuple(list_smiles)

    tup_values = tuple(
        "{0:.{1}f}\n".format(
            func(smiles, use_estimated_params=args.use_estimated_params),
            args.decimal,
        )
        for smiles in tup_smiles
    )

    if filepath_output:
        with open(filepath_output, mode="w", encoding="utf-8") as f:
            f.writelines(tup_values)
    else:
        sys.stdout.writelines(tup_values)


if __name__ == "__main__":
    main(sys.argv[1:], prog="askadskii")
