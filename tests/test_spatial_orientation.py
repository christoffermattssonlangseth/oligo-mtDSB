"""Test utils.spatial_orientation: rotated / mirrored copies of a section are brought back to the reference.

Run with pytest, or directly:  python tests/test_spatial_orientation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils import spatial_orientation as so  # noqa: E402

LAYOUT = {"Cortex": (0, 6), "Thalamus": (2, 2), "Hypothalamus": (2, -2), "Striatum": (-2, 1), "Hippocampus": (1, 4),
          "Fiber tracts": (0, 3), "Olfactory areas": (-3, -1), "Ventricular system": (2.5, 0), "Pallidum": (-1, -1)}


def section(sid, angle_deg, mirror, rng):
    rows = []
    a = np.radians(angle_deg)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    for reg, (x, y) in LAYOUT.items():
        p = rng.normal([x, y], 0.3, size=(200, 2))
        if mirror:
            p[:, 0] = -p[:, 0]
        p = p @ R.T * 100 + [5000, 3000]
        rows.append(pd.DataFrame({"sample_id": sid, "RBD_compartment_simplified": reg,
                                  "x_centroid": p[:, 0], "y_centroid": -p[:, 1]}))  # image y points down
    return pd.concat(rows)


def test_rotated_and_mirrored_sections_align():
    rng = np.random.default_rng(0)
    cells = pd.concat([section("ref", 0, False, rng), section("rot", 70, False, rng), section("mir", -120, True, rng)],
                      ignore_index=True)
    o = so.fit_orientation(cells, reference="ref")
    assert not o.loc["rot", "mirrored"] and o.loc["mir", "mirrored"]
    assert (o.rmsd < 20).all(), o.rmsd
    c = so.apply_orientation(cells, o)
    cen = c.groupby(["sample_id", "RBD_compartment_simplified"])[["x_oriented", "y_oriented"]].mean()
    for sid in ("rot", "mir"):
        assert np.allclose(cen.loc[sid].to_numpy(), cen.loc["ref"].to_numpy(), atol=15)
    ref = cen.loc["ref"]
    # canonical frame: thalamus straight above hypothalamus, cortex on top, hypothalamus right of striatum
    assert abs(ref.loc["Thalamus", "x_oriented"] - ref.loc["Hypothalamus", "x_oriented"]) < 15
    assert ref.loc["Cortex", "y_oriented"] > ref.loc["Hypothalamus", "y_oriented"]
    assert ref.loc["Hypothalamus", "x_oriented"] > ref.loc["Striatum", "x_oriented"]


if __name__ == "__main__":
    test_rotated_and_mirrored_sections_align()
    print("test_rotated_and_mirrored_sections_align: ok")
