"""Strict XYZ and explicit chemical metadata; atom indices are zero based."""
import json
from pathlib import Path
import numpy as np
from pyscf.data.elements import charge


def integer(value, name):
    if type(value) is not int:
        raise ValueError(f"{name} must be an integer")
    return value


def electron_check(symbols, q, spin):
    n = sum(charge(s) for s in symbols) - integer(q, "charge")
    integer(spin, "spin")
    if n < 0 or spin < 0 or spin > n or (n - spin) % 2:
        raise ValueError("electron count incompatible with spin=Nalpha-Nbeta")
    return n


def load(xyz, metadata):
    lines = Path(xyz).read_text().splitlines()
    n = int(lines[0])
    if n < 1 or len(lines) != n + 2:
        raise ValueError("XYZ must contain exactly the declared atom count and comment")
    symbols, coords = [], []
    for row in lines[2:]:
        parts = row.split()
        if len(parts) != 4 or charge(parts[0]) <= 0:
            raise ValueError("invalid XYZ atom row")
        symbols.append(parts[0])
        coords.append([float(v) for v in parts[1:]])
    if not np.isfinite(coords).all():
        raise ValueError("coordinates must be finite")
    m = json.loads(Path(metadata).read_text())
    electron_check(symbols, m["charge"], m["spin"])
    if not m.get("identity") or not m.get("protonation"):
        raise ValueError("identity and protonation metadata are required")
    seen = []
    fragments = m["fragments"]
    if not fragments:
        raise ValueError("at least one fragment required")
    for f in fragments:
        atoms = f["atoms"]
        if not atoms or any(type(i) is not int or i < 0 or i >= n for i in atoms):
            raise ValueError("invalid fragment atom indices")
        electron_check([symbols[i] for i in atoms], f["charge"], f["spin"])
        seen.extend(atoms)
    if sorted(seen) != list(range(n)):
        raise ValueError("each atom must belong to exactly one fragment")
    if sum(f["charge"] for f in fragments) != m["charge"]:
        raise ValueError("fragment charges must sum to complex charge")
    frozen = m.get("frozen_atoms", [])
    if len(set(frozen)) != len(frozen) or any(type(i) is not int or i < 0 or i >= n for i in frozen):
        raise ValueError("invalid frozen atom indices")
    for bond in m.get("bonds", []):
        if len(bond) != 3 or any(type(i) is not int or i < 0 or i >= n for i in bond[:2]) or bond[0] == bond[1]:
            raise ValueError("bonds must be [atom0, atom1, order]")
    return symbols, np.asarray(coords), m


def xyz_text(symbols, coords, comment="EDA_GPU; angstrom"):
    return str(len(symbols)) + "\n" + comment + "\n" + "".join(
        f"{s} {r[0]:.12f} {r[1]:.12f} {r[2]:.12f}\n" for s, r in zip(symbols, coords))
