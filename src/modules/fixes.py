"""The fixes of the reference's haplotype bugs that `prohap.py -fix_bugs 1` (config `fix_bugs: true`) applies
(IMPLEMENTATION_NOTES section 9). With -fix_bugs 0 none applies and the output is the reference's, byte for byte.

    from modules.fixes import selected
    fixes = selected(args.fix_bugs)      # a frozenset of names from ALL
    if 'conflicts' in fixes: ...

PROHAP_FIXES=name,name in the environment applies only the fixes it names, for development: tests of one fix at a
time, and the attribution of output differences to fixes (tests/fix_diff_report.py). It has no effect with
-fix_bugs 0.
"""
import os

# in the order they were added:
#   conflicts  variants on one haplotype copy conflict only when the bases they change overlap (section 9.2)
#   liftover   the liftover check compares alleles without their shared bases (an indel's anchor) (section 9.3)
#   stop       the change that creates the first stop codon is in the coding sequence (section 9.4)
#   start      a lost start codon is found again downstream instead of skipping the haplotype (section 9.5)
#   frameshift a frameshift is never synonymous (section 9.6)
ALL = ('conflicts', 'liftover', 'stop', 'start', 'frameshift')


def selected(fix_bugs, environ=None):
    """The fixes to apply: every one of ALL with fix_bugs on (or those in PROHAP_FIXES), none with it off."""
    if not fix_bugs:
        return frozenset()
    names = (os.environ if environ is None else environ).get('PROHAP_FIXES')
    if names is None:
        return frozenset(ALL)
    chosen = frozenset(n.strip() for n in names.split(',') if n.strip())
    unknown = sorted(chosen - set(ALL))
    if unknown:
        raise ValueError('PROHAP_FIXES names unknown fixes %s (known: %s)' % (', '.join(unknown), ', '.join(ALL)))
    return chosen
