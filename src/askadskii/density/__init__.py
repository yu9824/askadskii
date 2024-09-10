"""estimate density

Usage
-----

bash
^^^^
.. code-block:: bash
    askadskii density "*CC*"

Python
^^^^^^
.. code-block:: python
    from askadskii.density import estimate_density

    estimate_density("*CC*")

References
----------

- Askadskiĭ, A. A. (2003). Computational materials science of polymers. Cambridge Int Science Publishing.
- Askadskiĭ, A. A. (1996). Physical properties of polymers (Vol. 2). CRC Press.
"""

from ._core import estimate_density, estimate_vdw_volume

__all__ = ("estimate_density", "estimate_vdw_volume")

import askadskii.utils.dataset._decrypt  # noqa: F401
