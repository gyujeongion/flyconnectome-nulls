"""Seeds of the registered grids.

Every pre-registered grid in the paper was run on seeds 0-9. Later extensions put seeds 10-19
under the same run paths (runs/ecofix_P*, runs/dose_P*), so any script that globs those paths
would silently turn a registered n = 10 analysis into an n = 20 one the day the extension
lands. Scripts that reproduce registered numbers pass their globs through `registered`;
analyses of the extension import `EXTENSION` explicitly instead.
"""
import re

REGISTERED = range(10)
EXTENSION = range(10, 20)
_SEED = re.compile(r"_s(\d+)(?:/|$)")


def seed_of(path):
    m = _SEED.search(path)
    return int(m.group(1)) if m else None


def registered(paths):
    """keep only paths belonging to the registered seeds 0-9"""
    return [p for p in paths if seed_of(p) in REGISTERED]
