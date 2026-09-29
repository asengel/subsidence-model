from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from reservoir_subsidence.compaction import calculate_compaction_volume
from reservoir_subsidence.grid import make_surface_grid, results_dataframe
from reservoir_subsidence.van_opstal import calculate_subsidence
from reservoir_subsidence.visualization import (
    plot_subsidence_map,
    plot_center_profile,
)


def main():
    data = pd.read_csv(ROOT / "data" / "synthetic_reservoir.csv")

    data["compaction_volume_m3"] = calculate_compaction_volume(
        data["source_cell_area_m2"],
        data["thickness_m"],
        data["compaction_coefficient_per_bar"],
        data["pressure_depletion_bar"],
    )
    data["basement_depth_m"] = (
        data["reservoir_depth_m"] * data["basement_depth_factor"]
    )

    x, y, xx, yy = make_surface_grid(
        xmin_m=-4500,
        xmax_m=6500,
        ymin_m=-4000,
        ymax_m=4500,
        spacing_m=200,
    )

    reservoir_results = {}

    for reservoir_id, group in data.groupby("reservoir_id"):
        reservoir_results[reservoir_id] = calculate_subsidence(
            xx,
            yy,
            group["x_m"].to_numpy(),
            group["y_m"].to_numpy(),
            group["reservoir_depth_m"].to_numpy(),
            group["basement_depth_m"].to_numpy(),
            group["compaction_volume_m3"].to_numpy(),
        )

    total = sum(reservoir_results.values())

    result = results_dataframe(xx, yy, reservoir_results)
    result.to_csv(ROOT / "data" / "synthetic_subsidence_output.csv", index=False)

    plot_subsidence_map(
        xx,
        yy,
        total,
        output_path=ROOT / "figures" / "total_subsidence_map.png",
    )
    plot_center_profile(
        x,
        total,
        output_path=ROOT / "figures" / "subsidence_center_profile.png",
    )

    print(f"Maximum predicted subsidence: {total.max()*1000:.2f} mm")


if __name__ == "__main__":
    main()
