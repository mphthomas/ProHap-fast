"""The reference modules (legacy/src/modules), for the fast code to reuse.

Reusing the frozen functions (variant validity, conflicting variants, coordinates, population frequencies, FASTA
reading, ...) is what keeps the biology of ProHap-fast identical to the reference's:

    from modules.legacy_import import legacy
    legacy.vcf_reader.check_vcf_line_validity(...)
    legacy.coordinates_toolbox.get_rna_position(...)

Both trees name their package `modules` (a namespace package, no __init__.py), and the reference modules import each
other as `modules.<name>`. They are therefore imported once, with `modules` in sys.modules pointing at the reference
directory only, and src's own `modules` entries are put back afterwards. The reference modules keep their names
(`modules.common`, ...) as __name__ but live in sys.modules only under `prohap_legacy.<name>`.
"""
import importlib
import sys
import types
from pathlib import Path

LEGACY_MODULES = Path(__file__).resolve().parents[2] / 'legacy' / 'src' / 'modules'
NAMES = ('common', 'coordinates_toolbox', 'vcf_reader', 'get_haplotypes', 'process_haplotypes')
ALIAS = 'prohap_legacy'


def _load():
    if ALIAS in sys.modules:
        return sys.modules[ALIAS]
    own = {k: v for k, v in sys.modules.items() if k == 'modules' or k.startswith('modules.')}
    for k in own:
        del sys.modules[k]
    package = types.ModuleType('modules')
    package.__path__ = [str(LEGACY_MODULES)]
    package.__package__ = 'modules'
    sys.modules['modules'] = package
    alias = types.ModuleType(ALIAS)
    alias.__path__ = []
    alias.__doc__ = 'legacy/src/modules, imported by src/modules/legacy_import.py'
    try:
        for name in NAMES:
            module = importlib.import_module('modules.' + name)
            assert Path(module.__file__).resolve().parent == LEGACY_MODULES.resolve(), module.__file__
            setattr(alias, name, module)
    finally:
        loaded = {k: v for k, v in sys.modules.items() if k == 'modules' or k.startswith('modules.')}
        for k in loaded:
            del sys.modules[k]
        sys.modules.update(own)
    sys.modules[ALIAS] = alias
    for name in NAMES:
        sys.modules[ALIAS + '.' + name] = getattr(alias, name)
    return alias


legacy = _load()
