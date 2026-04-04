"""
Concrete class for running the ros functional test for FATES.
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from framework.functional_test import FunctionalTest
from framework.utils.plotting import blank_plot, get_color_palette

class FireMortTest(FunctionalTest):
    """Fire mortality test class"""
    name = 'fire_mortality'

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plot output associated with fuel tests

        Args:
            run_dir (str): run directory
            out_file (str): output file
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory
        """
        
        # get output file
        mortality_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))
        
        # observational datasets 
        obs_dat, sh_obs, tau_obs = self.get_obs_dfs()
        
        # plot observational comparisons
        self.bar_plot(mortality_dat, obs_dat, sh_obs, 'SH', 'fire_mortality_bySH', 
                      'Scorch Height (m)', save_figs, plot_dir)
        
        self.bar_plot(mortality_dat, obs_dat, tau_obs, 'tau_l', 'fire_mortality_bytau', 
                      'Fire Residence Time (min)', save_figs, plot_dir)
        
        self.plot_tau_c(mortality_dat, obs_dat, save_figs, plot_dir)
        
    @staticmethod
    def plot_tau_c(ds, obs_dat, save_figs, plot_dir):
        
        data_frame = pd.DataFrame({
            "dbh": np.tile(ds.dbh, len(ds.pft)),
            "pft": np.repeat(ds.pft, len(ds.dbh)),
            'tau_c': ds.tau_c.values.flatten()})
        
        max_dbh = data_frame["dbh"].max()
        max_tau = data_frame["tau_c"].max()
        
        blank_plot(max_dbh, 0.0, max_tau, 0.0, draw_horizontal_lines=True)
        
        pfts = np.unique(data_frame.pft.values)
        colors = get_color_palette(len(pfts))
        for rank, pft in enumerate(pfts):
            dat = data_frame[data_frame.pft == pft]
            plt.plot(
                dat.dbh.values,
                dat['tau_c'].values,
                lw=2,
                color=colors[rank],
                label=pft,
            )
        plt.scatter(obs_dat.diameter, obs_dat.tau_c, c='black', label='observations')
        plt.xlabel("DBH (cm)", fontsize=11)
        plt.ylabel("Critical Fire Residence Time for Cambial Burning (min)", fontsize=11)
        plt.title("Critical fire resdience time", fontsize=11)
        plt.legend(loc="upper left", title="PFT")
        if save_figs:
            fig_name = os.path.join(plot_dir, f"tauc_plot.png")
            plt.savefig(fig_name)
        
    @staticmethod
    def bar_plot(ds: xr.Dataset, tree_obs: pd.DataFrame, mortality_obs: pd.DataFrame, 
                 xvar: str, mortality_var: str, xlab: str, save_figs: bool, plot_dir: str):
        """Plots a bar plot based on input dataframe

        Args:
            df (xr.Dataset): input dataset
            tree_obs (pd.DataFrame): dataframe with observations about trees
            mortality_obs (pd.DataFrame): dataframe with observations about tree mortality
            xvar (str): variable to plot on x axis
            mortality_var (str): mortality variable
            xlab (str): xlabel
            save_figs (bool): whether or not to save figures
            plot_dir (str): where to save figures
        """
        
        data_frame = pd.DataFrame({
            "treeID": np.repeat(ds.treeID, len(ds[xvar])),
            xvar: np.tile(ds[xvar], len(ds.treeID)),
            'fire_mortality': ds[mortality_var].values.flatten()})
        data_frame['type'] = 'modeled'
        df = pd.concat([data_frame, mortality_obs])
        df = df.merge(tree_obs, on='treeID', how='inner')
        
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 6), sharey=True)
        axes = axes.flatten()

        bar_width = 0.3

        # unique values for x-axis
        x_values = np.unique(df[xvar])
        x = np.arange(len(x_values))

        legend_handles = []

        # loop over unique treeIDs and create subplots
        for i, treeID in enumerate(np.unique(df.treeID)):
            
            # extract data for this tree
            dat = df[df.treeID == treeID]  

            # extract modeled and observed fire mortality
            modeled = dat[dat['type'] == 'modeled']
            observed = dat[dat['type'] == 'observations']

            # plot bars side-by-side
            bars1 = axes[i].bar(x - bar_width/2, (modeled.fire_mortality)*100.0, 
                                width=bar_width, label='Modeled', color='blue')
            bars2 = axes[i].bar(x + bar_width/2, (observed.fire_mortality)*100.0, 
                                width=bar_width, label='Observed', color='orange')
            if i == 0:
                legend_handles.extend([bars1[0], bars2[0]])

            # for printing
            tree_dbh = observed.diameter.values[0]
            crown_length = observed.crown_length.values[0]
            height = observed.height.values[0]
            bark_thickness = observed.bark_thickness.values[0]
            
            # formatting
            axes[i].set_title(f"dbh: {tree_dbh} cm\nheight: {height} m\ncrown length: {crown_length} m\nbark thickness: {bark_thickness} cm")
            axes[i].set_xticks(x)
            axes[i].set_xticklabels(x_values) 
            axes[i].set_xlabel(xlab)
            axes[i].set_ylabel('Fire Mortality')

        fig.legend(legend_handles, ['Modeled', 'Observed'], loc='upper center', 
                   bbox_to_anchor=(0.5, 1.05), ncol=2)
        plt.tight_layout()
        if save_figs:
            fig_name = os.path.join(plot_dir, f"{xvar}_plot.png")
            plt.savefig(fig_name)

    
    @staticmethod
    def get_obs_dfs():
        """Return some hard-coded observational datasets
        Peterson & Ryan 1986 Environmental Management: 10(6)

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: output plots
        """
        
        obs_dat = pd.DataFrame({'species': ['douglas-fir', 'douglas-fir',
                                    'grand fir', 'grand fir',
                                    'subalpine fir', 'subalpine fir'],
                        'treeID': [1, 2, 3, 4, 5, 6],
                        'diameter': [20.0, 40.0, 20.0, 40.0, 20.0, 40.0],
                        'height': [15.5, 24.4, 17.4, 30.7, 15.1, 24.1],
                        'crown_length': [12.0, 18.5, 13.5, 23.1, 11.9, 18.2],
                        'bark_thickness': [1.1, 2.2, 0.9, 1.7, 0.3, 0.6],
                        'tau_c': [3.2, 13.4, 2.1, 8.3, 0.3, 1.0]})
        
        sh_obs = pd.DataFrame({'fire_mortality': [0.94, 0.0, 1.0, 0.0, 1.0, 1.0,
                                          0.99, 0.03, 1.0, 0.13, 1.0, 1.0,
                                          1.0, 0.72, 1.0, 0.71, 1.0, 1.0],
                      'SH': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                             10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                             20.0, 20.0, 20.0, 20.0, 20.0, 20.0],
                      'treeID': [1, 2, 3, 4, 5, 6,
                                1, 2, 3, 4, 5, 6,
                                1, 2, 3, 4, 5, 6]})
        sh_obs['type'] = 'observations'
        
        tau_obs = pd.DataFrame({'fire_mortality': [0.74, 0.0, 0.82, 0.0, 1.0, 1.0,
                                          0.83, 0.0, 0.9, 0.01, 1.0, 1.0,
                                          0.99, 0.03, 1.0, 0.13, 1.0, 1.0],
                      'tau_l': [2.4, 2.4, 2.4, 2.4, 2.4, 2.4,
                             3.1, 3.1, 3.1, 3.1, 3.1, 3.1,
                             6, 6, 6, 6, 6, 6],
                      'treeID': [1, 2, 3, 4, 5, 6,
                                1, 2, 3, 4, 5, 6,
                                1, 2, 3, 4, 5, 6]})
        tau_obs['type'] = 'observations'
        
        return obs_dat, sh_obs, tau_obs
        
    
