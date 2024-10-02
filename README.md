# askadskii

Askadskii's atomic group contribution method.

## Usage

```plaintext
usage: askadskii density [-h] [-i INPUT_FILE] [-d DECIMAL] [--debug] [--deprecated] [-o OUTPUT_FILE] [-e] [SMILES ...]

positional arguments:
  SMILES                smiles

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        A filepath containing a list of SMILES separated by newlines
  -d DECIMAL, --decimal DECIMAL
                        the number of decimal places to display, by default 3
  --debug               debug mode
  --deprecated          use deprecated method
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        filepath output
  -e, --use-estimated-params
                        use estimated covalent bond dinstance
```

## Install

You can install with `pip` ;

```bash
pip install git+https://github.com/yu9824/askadskii.git

```

## References

- Askadskiĭ, A. A. (2003). Computational materials science of polymers. Cambridge Int Science Publishing.
- Askadskiĭ, A. A. (1996). Physical properties of polymers (Vol. 2). CRC Press.
