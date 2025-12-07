
"""
Concrete class for running the canopyfuel functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from functional_class import FunctionalTest
from utils_plotting import get_color_palette

class CanopyFuelTest(FunctionalTest):
    """Canopy fuel test class"""

    name = "canopyfuel"

    def __init__(self, test_dict):
        super().__init__(
            CanopyFuelTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plot output associated with canopy fuel tests

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """

        # read in canopy fuel data
        cfuel_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        ## 1. some barchart plots for variables with one dimension

        self.plot_barchart(
            cfuel_dat,
            "CBD",
            "Canopy bulk density",
            "kg m$^{-3}$",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "CBH",
            "Canopy base height",
            "m",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
        self.plot_barchart(
            cfuel_dat,
            "ROS_min",
            "Critical rate of spread",
            "m min$^{-1}$",
            save_figs,
            plot_dir,
            by_fuel_model=False,
            stacked= False,
        )
    
        ## 2. plots for variables with at least 2 dimensions
        self.plot_facet(
            cfuel_dat, 
            x="canopy_water", y="FI_init",
            hue="patch_type",row=None,col=None,style=None,
            kind="line", agg=None,
            x_label="Canopy water content(%)",
            y_label="Fire intensity to ignite crown fuel(kW m-1)",
            save_figs=save_figs,
            plot_dir=plot_dir,
            filename="FI_init_plot.png",
            )
        self.plot_facet(
           cfuel_dat,
           x="wind", y="ROS_front",
           hue="fire_weather",row="fuel_model",col="patch_type",style=None,
           kind="line", agg=None,
           x_label="Wind speed(m min-1)",
           y_label="Surface fire rate of spread(m min-1)",
           save_figs=save_figs,
           plot_dir=plot_dir,
           filename="ROS_front_plot.png",
        )
        self.plot_facet(
            cfuel_dat,
            x="wind", y="ROS_active",
            hue="fire_weather",row="canopy_water",col="patch_type",style="fuel_model",
            kind="line", agg=None,
            x_label="Wind speed(m min-1)",
            y_label="Active crown fire rate of spread(m min-1)",
            save_figs=save_figs,
            plot_dir=plot_dir,
            filename="ROS_act_plot.png",
        )
        self.plot_facet(
            cfuel_dat,
            x="wind",y="ROS_final",
            hue="fire_weather",row="canopy_water",col="patch_type",style="fuel_model",
            kind="line", agg=None,
            x_label="Wind speed(m min-1)",
            y_label="Final rate of spread(m min-1)",
            save_figs=save_figs,
            plot_dir=plot_dir,
            filename="ROS_final_plot.png",
        )
        self.plot_facet(
            cfuel_dat,
            x="wind", y="FI_final",
            hue="fire_weather",row="canopy_water",col="patch_type",style="fuel_model",
            kind="line", agg=None,
            x_label="Wind speed(m min-1)",
            y_label="Final fire intensity(kW m-1)",
            save_figs=save_figs,
            plot_dir=plot_dir,
            filename="FI_final_plot.png",
        )
        self.plot_facet(
            cfuel_dat,
            x="wind",y="CFB",
            hue="fire_weather",row="canopy_water",col="patch_type",style="fuel_model",
            kind="line", agg=None,
            x_label="Wind speed(m min-1)",
            y_label="Crown fraction burnt",
            save_figs=save_figs,
            plot_dir=plot_dir,
            filename="CFB_plot.png",

        )

# ------- Define plot functions -------

    @staticmethod
    def plot_barchart(
        cfuel_dat: xr.Dataset,
        var: str,
        varname: str,
        units: str,
        save_figs: bool,
        plot_dir: str,
        by_fuel_model: bool = False,
        stacked: bool = False,
          ):
        """ Plot canopy fuel test outputs as bar plot
        Args:
        cfuel_dat (xr.Dataset): canopy fuel test data output
        var (str): variable to plot
        varname (str): variable name for y-axis
        units (str): unit expression
        save_figs (str): dir to save figure
        plot_dir (str): which dir to save figs
        by_fuel_model (bool, optional): whether or not the bar plot is grouped by fuel model, default to True

        """
        assert var in cfuel_dat, f"{var!r} not found in dataset variables."
        da = cfuel_dat[var]

        # figure out the dimensions we care about
        has_patch = "patch_type" in da.dims
        has_fm    = "fuel_model" in da.dims

        if not has_patch:
            raise ValueError(f"{var} is missing 'patch_type' dimension; got dims {da.dims}")
    
        patch_types = [str(p) for p in cfuel_dat["patch_type"].values]


        if by_fuel_model and has_fm:
            da2 = da.transpose("patch_type", "fuel_model")
            arr = da2.to_numpy()
            fm_labels = [str(f) for f in cfuel_dat["fuel_model"].values]

            fig, ax = plt.subplots(figsize=(8, 4.6))
            x = np.arange(len(patch_types))
            n_fm = arr.shape[1]
            width = 0.8 / n_fm

            if stacked:
                bottom = np.zeros(len(patch_types))
                for j in range(n_fm):
                    ax.bar(
                        x,
                        arr[:, j],
                        width=0.8, 
                        bottom=bottom,
                        label=fm_labels[j],
                    )
                    bottom += arr[:, j]
            else:
                for j in range(n_fm):
                    offset = (j - (n_fm - 1) / 2) * width
                    ax.bar(
                        x + offset,
                        arr[:, j],
                        width=width,
                        label=fm_labels[j],
                        )
            ax.set_xticks(x)
            ax.set_xticklabels(patch_types, rotation=90)
            ax.legend(title="fuel_model", loc="center left", bbox_to_anchor=(1, 0.5))
        else:
            if has_fm:
                da_plot = da.mean(dim="fuel_model")
            else:
                da_plot = da
            
            da_plot = da_plot.transpose("patch_type", ...)
            y = da_plot.to_numpy()
            fig, ax = plt.subplots(figsize=(7.5, 4.2))
            ax.bar(patch_types, y)
            ax.set_xticklabels(patch_types, rotation=90)

        ax.set_ylabel(f"{varname} ({units})", fontsize=11)
        ax.set_xlabel("Patch Type")
        ax.set_title(varname)
        plt.tight_layout()

        if save_figs:
            fig_name = os.path.join(plot_dir, f"{varname}_plot.png")
            plt.savefig(fig_name)

    @staticmethod
    def plot_facet(
        ds, *,
        x: str,
        y: str,
        hue: str | None = None,
        row: str | None = None,
        col: str | None = None,
        style: str | None = None,
        kind: str = "line",      # "line" or "scatter"
        agg: str | None = "mean",# "mean", "median", or None (no aggregation)
        orders: dict | None = None,  # e.g., {"CWC":[40,50,60,70,80], "NI":[1000,2000,...]}
        x_label: str | None = None,
        y_label: str | None = None,
        title: str | None = None,
        save_figs: bool = False,
        plot_dir: str | None = None,
        filename: str = "facet_plot.png",
    ):
        """
        General faceted plot using matplotlib.

        Panels: row × col
        Inside each panel: one series per `hue`; optionally split each hue by `style`.

        - If `agg` is provided, aggregates y vs x within each group (good for many points).
        - If `kind="scatter"`, plots points; if "line", draws lines (with optional markers).
        """

        # Tidy DataFrame
        if hasattr(ds, "to_dataframe"):
            df = ds.to_dataframe().reset_index()
        else:
            df = ds.copy()

        # Collect the fields actually used
        fields = [x, y] + [v for v in [hue, row, col, style] if v]
        missing = [k for k in fields if k not in df.columns]
        if missing:
            raise KeyError(f"Missing columns: {missing}")

        df = df.dropna(subset=[x, y])  # basic hygiene

        # Helpers for level orders
        orders = orders or {}
        def _levels(key):
            return orders.get(key, sorted(pd.unique(df[key]))) if key else [None]

        row_levels = _levels(row)
        col_levels = _levels(col)
        hue_levels = _levels(hue)

        n_rows, n_cols = len(row_levels), len(col_levels)
        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(3.8 * n_cols, 2.8 * n_rows),
            sharex=True, sharey=True
        )
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = np.array([axes])
        elif n_cols == 1:
            axes = axes[:, None]

        # Aggregator
        if agg is None:
            def _agg(g):
                return g.sort_values(x)[[x, y]]
        else:
            agg_name = {"mean": "mean", "median": "median"}[agg]
            def _agg(g):
                out = g.groupby(x, as_index=False)[y].agg(agg_name)
                return out.sort_values(x)

        # Plot panels
        for r, rv in enumerate(row_levels):
            for c, cv in enumerate(col_levels):
                ax = axes[r, c]
                sub = df.copy()
                if row: sub = sub[sub[row] == rv]
                if col: sub = sub[sub[col] == cv]

                ttl_bits = []
                if row: ttl_bits.append(f"{row}={rv}")
                if col: ttl_bits.append(f"{col}={cv}")
                ax.set_title(" • ".join(ttl_bits) if ttl_bits else (title or ""))

                if sub.empty:
                    ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
                    continue

                # Loop hues
                for hv in hue_levels:
                    sub_h = sub if not hue else sub[sub[hue] == hv]
                    if sub_h.empty:
                        continue

                    # Optional style split
                    if style:
                        for sv, sub_s in sub_h.groupby(style):
                            cur = _agg(sub_s)
                            if cur.empty: 
                                continue
                            if kind == "scatter":
                                ax.scatter(cur[x].to_numpy(), cur[y].to_numpy(),
                                        label=f"{hue}={hv}, {style}={sv}" if hue else f"{style}={sv}", s=18)
                            else:
                                ax.plot(cur[x].to_numpy(), cur[y].to_numpy(),
                                        marker="o", linestyle="-", linewidth=1.25, markersize=3.5,
                                        label=f"{hue}={hv}, {style}={sv}" if hue else f"{style}={sv}")
                    else:
                        cur = _agg(sub_h)
                        if cur.empty:
                            continue
                        if kind == "scatter":
                            ax.scatter(cur[x].to_numpy(), cur[y].to_numpy(),
                                    label=f"{hue}={hv}" if hue else None, s=18)
                        else:
                            ax.plot(cur[x].to_numpy(), cur[y].to_numpy(),
                                    marker="o", linestyle="-", linewidth=1.25, markersize=3.5,
                                    label=f"{hue}={hv}" if hue else None)

                if r == n_rows - 1:
                    ax.set_xlabel(x_label or x)
                if c == 0:
                    ax.set_ylabel(y_label or y)

        # Consolidated legend
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5),
                    title=(" × ".join([k for k in [hue, style] if k])) or None)

        if title and (n_rows * n_cols > 1):
            fig.suptitle(title, y=1.02)

        fig.tight_layout(rect=(0, 0, 0.9, 1))
        if save_figs and plot_dir:
            os.makedirs(plot_dir, exist_ok=True)
            fig.savefig(os.path.join(plot_dir, filename), dpi=150, bbox_inches="tight")
        plt.close(fig)

  
