"""
demo.py - simplified illustration of 3D bin packing with rotations,
weight limit and stacking support.

IMPORTANT
---------
This file is NOT the original solver. It is a small, self-contained
re-implementation written only to illustrate the problem and the basic
construction idea (greedy placement on corner points with rotations).
The full multi-stage heuristic (instance analysis, multi-ordering
construction, local search, exact sub-models) is kept private, see README.md.

What this demo does
-------------------
1. Generates a synthetic set of rectangular items with a fixed random seed.
2. Packs them one by one into identical bins using a greedy rule:
   lowest-then-frontmost-then-leftmost free corner, first allowed
   orientation that fits, respecting bounds, no overlap, weight limit and a
   minimum supported fraction of the item's base.
3. Prints the number of bins used and the volume utilisation of each bin.
4. Saves 3D renderings of the first two bins to images/.

Dependencies: numpy, matplotlib only.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

# ---------------------------------------------------------------------------
# Parameters of the synthetic instance
# ---------------------------------------------------------------------------
SEED = 7
N_ITEMS = 60
BIN_SIZE = (120.0, 80.0, 100.0)   # (x, y, z) inner dimensions, arbitrary units
BIN_MAX_WEIGHT = 900.0
MIN_SUPPORT = 0.70                # fraction of base area that must rest on something
EPS = 1e-9

# Item "families": (dimensions, weight range, allowed orientations)
# Orientation codes are the six permutations of (w, d, h), see ORIENTATIONS.
ITEM_FAMILIES = [
    # name,      (w, d, h),        weight,    allowed orientations
    ("crate",    (40.0, 30.0, 25.0), (20, 40), (0, 1)),            # upright only
    ("box",      (30.0, 20.0, 20.0), (8, 18),  (0, 1, 2, 3, 4, 5)),  # free
    ("slab",     (60.0, 40.0, 10.0), (15, 30), (0, 1)),            # flat only
    ("tube",     (20.0, 20.0, 60.0), (10, 25), (0, 1, 2, 3, 4, 5)),  # free
    ("cube",     (25.0, 25.0, 25.0), (12, 22), (0,)),              # fixed
]

# Colours per family (fixed order, one hue per family).
FAMILY_COLOURS = {
    "crate": "#2a78d6",
    "box":   "#eb6834",
    "slab":  "#1baf7a",
    "tube":  "#eda100",
    "cube":  "#e87ba4",
}

# The six axis-aligned orientations as permutations of (w, d, h).
# ORIENTATIONS[k] = indices of (w, d, h) placed along (x, y, z).
ORIENTATIONS = [
    (0, 1, 2),  # 0: w d h
    (1, 0, 2),  # 1: d w h   (rotation about z)
    (0, 2, 1),  # 2: w h d   (rotation about x)
    (2, 1, 0),  # 3: h d w   (rotation about y)
    (1, 2, 0),  # 4: d h w
    (2, 0, 1),  # 5: h w d
]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class Item:
    id: int
    family: str
    dims: tuple[float, float, float]   # nominal (w, d, h)
    weight: float
    allowed: tuple[int, ...]

    @property
    def volume(self) -> float:
        w, d, h = self.dims
        return w * d * h

    def extents(self, orient: int) -> tuple[float, float, float]:
        """Dimensions along (x, y, z) after applying an orientation code."""
        perm = ORIENTATIONS[orient]
        return tuple(self.dims[p] for p in perm)  # type: ignore[return-value]


@dataclass
class Placement:
    item: Item
    pos: tuple[float, float, float]
    ext: tuple[float, float, float]
    orient: int


@dataclass
class Bin:
    size: tuple[float, float, float]
    max_weight: float
    placements: list[Placement] = field(default_factory=list)
    corners: list[tuple[float, float, float]] = field(
        default_factory=lambda: [(0.0, 0.0, 0.0)]
    )
    weight: float = 0.0

    @property
    def volume(self) -> float:
        return self.size[0] * self.size[1] * self.size[2]

    @property
    def used_volume(self) -> float:
        return sum(p.item.volume for p in self.placements)

    @property
    def utilisation(self) -> float:
        return self.used_volume / self.volume


# ---------------------------------------------------------------------------
# Feasibility checks (mirror the constraints described in the README)
# ---------------------------------------------------------------------------
def inside(bin_: Bin, pos, ext) -> bool:
    return all(pos[i] + ext[i] <= bin_.size[i] + EPS for i in range(3))


def overlaps(a_pos, a_ext, b_pos, b_ext) -> bool:
    for i in range(3):
        if a_pos[i] + a_ext[i] <= b_pos[i] + EPS:
            return False
        if b_pos[i] + b_ext[i] <= a_pos[i] + EPS:
            return False
    return True


def supported_fraction(bin_: Bin, pos, ext) -> float:
    """Fraction of the base area resting on the floor or on item tops."""
    if pos[2] <= EPS:
        return 1.0
    base_area = ext[0] * ext[1]
    covered = 0.0
    for p in bin_.placements:
        top = p.pos[2] + p.ext[2]
        if abs(top - pos[2]) > 1e-6:
            continue
        dx = min(pos[0] + ext[0], p.pos[0] + p.ext[0]) - max(pos[0], p.pos[0])
        dy = min(pos[1] + ext[1], p.pos[1] + p.ext[1]) - max(pos[1], p.pos[1])
        if dx > 0 and dy > 0:
            covered += dx * dy
    return covered / base_area


def feasible(bin_: Bin, item: Item, pos, ext) -> bool:
    if bin_.weight + item.weight > bin_.max_weight + EPS:
        return False
    if not inside(bin_, pos, ext):
        return False
    for p in bin_.placements:
        if overlaps(pos, ext, p.pos, p.ext):
            return False
    return supported_fraction(bin_, pos, ext) >= MIN_SUPPORT - 1e-9


# ---------------------------------------------------------------------------
# Greedy construction
# ---------------------------------------------------------------------------
def try_place(bin_: Bin, item: Item) -> bool:
    """Scan candidate corners from the bottom, front, left; first fit wins."""
    for corner in sorted(bin_.corners, key=lambda c: (c[2], c[1], c[0])):
        for orient in item.allowed:
            ext = item.extents(orient)
            if feasible(bin_, item, corner, ext):
                bin_.placements.append(Placement(item, corner, ext, orient))
                bin_.weight += item.weight
                x, y, z = corner
                for new in ((x + ext[0], y, z), (x, y + ext[1], z), (x, y, z + ext[2])):
                    if new not in bin_.corners and inside(bin_, new, (0, 0, 0)):
                        bin_.corners.append(new)
                return True
    return False


def pack(items: list[Item]) -> list[Bin]:
    """First-fit-decreasing by volume over identical bins."""
    bins: list[Bin] = []
    for item in sorted(items, key=lambda it: -it.volume):
        if not any(try_place(b, item) for b in bins):
            new_bin = Bin(BIN_SIZE, BIN_MAX_WEIGHT)
            if not try_place(new_bin, item):
                raise RuntimeError(f"item {item.id} does not fit in an empty bin")
            bins.append(new_bin)
    return bins


def verify(bins: list[Bin], items: list[Item]) -> None:
    """Independent re-check of the final layout against every constraint."""
    placed_ids = [p.item.id for b in bins for p in b.placements]
    assert sorted(placed_ids) == sorted(it.id for it in items), "items missing or duplicated"
    for b in bins:
        assert sum(p.item.weight for p in b.placements) <= b.max_weight + EPS, "weight limit"
        for i, p in enumerate(b.placements):
            assert p.orient in p.item.allowed, "forbidden orientation"
            assert p.ext == p.item.extents(p.orient), "extent mismatch"
            assert all(p.pos[k] >= -EPS for k in range(3)) and inside(b, p.pos, p.ext), "bounds"
            for q in b.placements[i + 1:]:
                assert not overlaps(p.pos, p.ext, q.pos, q.ext), "overlap"
            others = Bin(b.size, b.max_weight, [q for q in b.placements if q is not p])
            assert supported_fraction(others, p.pos, p.ext) >= MIN_SUPPORT - 1e-9, "support"


# ---------------------------------------------------------------------------
# Synthetic instance
# ---------------------------------------------------------------------------
def generate_items(rng: np.random.Generator, n: int) -> list[Item]:
    items = []
    for i in range(n):
        name, dims, (w_lo, w_hi), allowed = ITEM_FAMILIES[rng.integers(len(ITEM_FAMILIES))]
        items.append(Item(i, name, dims, float(rng.uniform(w_lo, w_hi)), allowed))
    return items


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def box_faces(pos, ext):
    x, y, z = pos
    dx, dy, dz = ext
    v = np.array([
        [x, y, z], [x + dx, y, z], [x + dx, y + dy, z], [x, y + dy, z],
        [x, y, z + dz], [x + dx, y, z + dz], [x + dx, y + dy, z + dz], [x, y + dy, z + dz],
    ])
    return [
        [v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]],
        [v[0], v[1], v[5], v[4]], [v[2], v[3], v[7], v[6]],
        [v[1], v[2], v[6], v[5]], [v[0], v[3], v[7], v[4]],
    ]


def render_bin(bin_: Bin, index: int, path: str) -> None:
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    for p in bin_.placements:
        faces = Poly3DCollection(
            box_faces(p.pos, p.ext),
            facecolors=FAMILY_COLOURS[p.item.family],
            edgecolors="#333333", linewidths=0.6, alpha=0.9,
        )
        ax.add_collection3d(faces)
    sx, sy, sz = bin_.size
    ax.set_xlim(0, sx); ax.set_ylim(0, sy); ax.set_zlim(0, sz)
    ax.set_box_aspect((sx, sy, sz))
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title(
        f"Bin {index}: {len(bin_.placements)} items, "
        f"volume utilisation {bin_.utilisation:.1%}, "
        f"weight {bin_.weight:.0f}/{bin_.max_weight:.0f}"
    )
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=c, label=f) for f, c in FAMILY_COLOURS.items()
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    ax.view_init(elev=25, azim=-55)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    rng = np.random.default_rng(SEED)
    items = generate_items(rng, N_ITEMS)
    bins = pack(items)
    verify(bins, items)

    total_item_volume = sum(it.volume for it in items)
    lower_bound = math.ceil(total_item_volume / bins[0].volume - 1e-9)

    print("Simplified 3D bin packing demo (not the original solver)")
    print(f"seed={SEED}  items={N_ITEMS}  bin={BIN_SIZE}  max weight={BIN_MAX_WEIGHT}  "
          f"min support={MIN_SUPPORT:.0%}")
    print(f"bins used: {len(bins)}   (volume lower bound: {lower_bound})")
    for i, b in enumerate(bins, start=1):
        print(f"  bin {i}: {len(b.placements):2d} items, "
              f"volume utilisation {b.utilisation:6.1%}, "
              f"weight {b.weight:6.1f}/{b.max_weight:.0f}")
    overall = total_item_volume / (len(bins) * bins[0].volume)
    print(f"overall volume utilisation: {overall:.1%}")
    print("self-check: all constraints verified (bounds, overlap, orientation, weight, support)")

    os.makedirs("images", exist_ok=True)
    for i, b in enumerate(bins[:2], start=1):
        path = os.path.join("images", f"bin_{i}.png")
        render_bin(b, i, path)
        print(f"saved {path}")


if __name__ == "__main__":
    main()
