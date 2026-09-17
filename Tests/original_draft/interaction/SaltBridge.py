"""
This module is part of DuIvyProcedures.procedures, designed for dealing Sltbr. 
Written by 杜艾维.
"""

import os
import sys
import numpy as np
import seaborn as sns
import MDAnalysis as mda
from itertools import chain
from scipy.integrate import simpson
from matplotlib import pyplot as plt
from MDAnalysis.analysis import distances as mda_dist
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM

base = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
if base not in sys.path:
    sys.path.insert(0, base)

from utils import log
from framework.confParser import Config


class SaltBridge(log):
    def __init__(self, config: Config) -> None:
        self.conf = config
        self.load_mplstyle()

    @log.module_decorator
    def __call__(self) -> None:

        tpr = self.conf["tpr"]
        xtc = self.conf["xtc"]

        dist_cutoff = self.conf["dist_cutoff"]
        byIndex = self.conf["byIndex"]
        group = self.conf["group"]
        ignore_chain_end = self.conf["ignore_chain_end"]
        positive_Index = self.conf["positive_Index"]
        negative_Index = self.conf["negative_Index"]
        NH3_atomnames = self.conf["NH3_atomnames"]
        COO_atomnames = self.conf["COO_atomnames"]
        Backbone_atomnames = self.conf[
            "Backbone_atomnames"
        ]  #  ["H", "N", "CA", "C", "O"]
        positive_Index = np.array(positive_Index) - 1  # mda atom index (start from 0)
        negative_Index = np.array(negative_Index) - 1
        calc_lifetime = self.conf["calc_lifetime"]
        tau_max = self.conf["tau_max"]
        window_step = self.conf["window_step"]
        intermittency = self.conf["intermittency"]

        fstart = self.conf.get("frame_start", None)
        fend = self.conf.get("frame_end", None)
        fstep = self.conf.get("frame_step", None)
        for key in [fstart, fend, fstep]:
            if key is None:
                continue
            if isinstance(key, int):
                if key < 0:
                    self.critical(
                        f"frame_start, frame_end, frame_step should not be negative integers, but got {key}"
                    )
            else:
                self.critical(
                    f"frame_start, frame_end, frame_step should be integers or leave blank, but got {key} of type {type(key)}"
                )
        self.info(f"Analyzing on trajectory[{fstart}:{fend}:{fstep}]")

        u = mda.Universe(f"../{tpr}", f"../{xtc}")
        if byIndex == True:
            all_index = [i for i in chain(*positive_Index, *negative_Index)]
            if len(all_index) != 0:
                all_atoms = u.atoms[all_index]
                all_atoms.write("SaltBridges_byIndex.pdb")
            else:
                self.error("Empty Indexs, please check your input.")
                return
            time_array, dist_dataframe = [], []
            for ts in u.trajectory[fstart:fend:fstep]:
                positive_centers, negative_centers = (
                    [],
                    [],
                )
                for AA_indexs in positive_Index:
                    atoms = u.atoms[AA_indexs]
                    atom_centers = atoms.center_of_mass()
                    positive_centers.append(atom_centers)
                for AA_indexs in negative_Index:
                    atoms = u.atoms[AA_indexs]
                    atom_centers = atoms.center_of_mass()
                    negative_centers.append(atom_centers)
                positive_centers = np.array(positive_centers)
                negative_centers = np.array(negative_centers)
                dist_matrix = mda_dist.distance_array(
                    positive_centers, negative_centers, box=ts.dimensions
                )
                dist_matrix *= 0.1  ## A to nm
                time_array.append(ts.time)
                dist_dataframe.append(dist_matrix)
            time_array = np.array(time_array)
            dist_dataframe = np.array(dist_dataframe)

            positive_residues, negative_residues = [], []
            for AA_indexs in positive_Index:
                res = u.atoms[AA_indexs].residues
                positive_residues.append(f"{res.resnames[0]}{res.resnums[0]}")
            for AA_indexs in negative_Index:
                res = u.atoms[AA_indexs].residues
                negative_residues.append(f"{res.resnames[0]}{res.resnums[0]}")

        else:  # by charge !
            positive_uAAs, negative_uAAs = [], []
            positive_residues, negative_residues = [], []
            positive_uAAc, negative_uAAc = [], []
            positive_residues_c, negative_residues_c = [], []
            for res in u.select_atoms(group).residues:
                res_set = set(res.atoms.names)
                iswhole = set(Backbone_atomnames).issubset(res_set)
                if iswhole == True:  ## to obtain charge center of sidechain
                    key = " and ".join([f"not name {a}" for a in Backbone_atomnames])
                    side = res.atoms.select_atoms(key)
                    if np.sum(side.charges) > 0.42:
                        positive_uAAs.append(side)
                        positive_residues.append(f"{res.resname}{res.resnum}")
                    elif np.sum(side.charges) < -0.42:
                        negative_uAAs.append(side)
                        negative_residues.append(f"{res.resname}{res.resnum}")
                else:  ## not whole, get C-ter, N-ter, and sidechains
                    if set(NH3_atomnames).issubset(res_set):
                        key = " or ".join([f"name {a}" for a in NH3_atomnames])
                        key += " or name CA"  # add CA, which always charged
                        NH3_atoms = res.atoms.select_atoms(key)
                        NH3_name = f"{res.resname}{res.resnum}_NH3"
                        if np.sum(NH3_atoms.charges) > 0.42:
                            positive_uAAc.append(NH3_atoms)
                            positive_residues_c.append(NH3_name)
                        else:
                            self.warn(
                                f"The NH3 group of {NH3_name} contained charges less than 0.42, quite wired ? DIP ignored it."
                            )
                            print(NH3_atoms.charges)
                    if set(COO_atomnames).issubset(res_set):
                        key = " or ".join([f"name {a}" for a in COO_atomnames])
                        key += " or name CA"  # add CA, which always charged
                        COO_atoms = res.atoms.select_atoms(key)
                        COO_name = f"{res.resname}{res.resnum}_COO"
                        if np.sum(COO_atoms.charges) < -0.42:
                            negative_uAAc.append(COO_atoms)
                            negative_residues_c.append(COO_name)
                        else:
                            self.warn(
                                f"The COO group of {COO_name} contained charges higher than -0.42, quite wired ? DIP ignored it."
                            )
                            print(COO_atoms.charges)
                    NH3_COO_Backbone = set(
                        NH3_atomnames + COO_atomnames + Backbone_atomnames
                    )
                    key = " and ".join([f"not name {a}" for a in NH3_COO_Backbone])
                    side = res.atoms.select_atoms(key)
                    if np.sum(side.charges) > 0.42:
                        positive_uAAc.append(side)
                        positive_residues_c.append(f"{res.resname}{res.resnum}_side")
                    elif np.sum(side.charges) < -0.42:
                        negative_uAAc.append(side)
                        negative_residues_c.append(f"{res.resname}{res.resnum}_side")

            self.info(
                f"Found {len(positive_uAAs)} AA with POSITIVE sidechain and {len(negative_uAAs)} AA with NEGATIVE sidechain in system"
            )
            self.info(
                f"And found {len(positive_uAAc)} AA (broken backbone) with POSITIVE charge and {len(negative_uAAc)} AA (broken backbone) with NEGATIVE charge in system"
            )
            if not ignore_chain_end:
                positive_uAAs += positive_uAAc
                negative_uAAs += negative_uAAc
                positive_residues += positive_residues_c
                negative_residues += negative_residues_c

            positive_Index = [side.ids for side in positive_uAAs]
            negative_Index = [side.ids for side in negative_uAAs]
            all_index = [i for i in chain(*positive_Index, *negative_Index)]
            if len(all_index) != 0:
                all_atoms = u.atoms[all_index]
                all_atoms.write("SaltBridges_byDefault.pdb")
            else:
                self.error("Empty Indexs, No charged groups found in your structure.")
                return

            time_array, dist_dataframe = [], []
            for ts in u.trajectory[fstart:fend:fstep]:
                positive_centers, negative_centers = (
                    [],
                    [],
                )
                for side in positive_uAAs:
                    atom_centers = side.center_of_charge()
                    positive_centers.append(atom_centers)
                for side in negative_uAAs:
                    atom_centers = side.center_of_charge()
                    negative_centers.append(atom_centers)
                positive_centers = np.array(positive_centers)
                negative_centers = np.array(negative_centers)
                dist_matrix = mda_dist.distance_array(
                    positive_centers, negative_centers, box=ts.dimensions
                )
                dist_matrix *= 0.1  ## A to nm
                time_array.append(ts.time)
                dist_dataframe.append(dist_matrix)
            time_array = np.array(time_array)
            dist_dataframe = np.array(dist_dataframe)

        ## output atoms indexs for checking and reuse this module
        with open("SaltBridge_Indexs.txt", "w") as fo:
            fo.write("PositiveGroups, Indexs\n")
            for i, index in enumerate(positive_Index):
                fo.write(f"{positive_residues[i]}, {[id+1 for id in index]}\n")
            fo.write("NegativeGroups, Indexs\n")
            for i, index in enumerate(negative_Index):
                fo.write(f"{negative_residues[i]}, {[id+1 for id in index]}\n")

        ## analyze results
        t_num, x_num, y_num = dist_dataframe.shape
        sltbr_datas, sltbr_issltbr = [], []
        sltbr_names, sltbr_occs, sltbr_dists = [], [], []
        sltbr_dist_ave_map = np.zeros((x_num, y_num))
        sltbr_dist_std_map = np.zeros((x_num, y_num))
        sltbr_occ_map = np.zeros((x_num, y_num))
        for x in range(x_num):
            for y in range(y_num):
                data = dist_dataframe[:, x, y]
                issltbr = data < dist_cutoff
                occupancy = np.sum(issltbr) / t_num
                sltbr_occ_map[x, y] = occupancy
                if occupancy > 0:
                    sltbr_datas.append(data)
                    sltbr_issltbr.append(issltbr)
                    sltbr_name = f"{positive_residues[x]}-{negative_residues[y]}"
                    sltbr_names.append(sltbr_name)
                    sltbr_occs.append(occupancy)
                    data_issltbr = data[issltbr]
                    dist_ave = np.average(data_issltbr)
                    dist_std = np.std(data_issltbr, ddof=1)
                    sltbr_dists.append((dist_ave, dist_std))
                    sltbr_dist_ave_map[x, y] = dist_ave
                    sltbr_dist_std_map[x, y] = dist_std

        # output info to csv
        with open("SaltBridge_info.csv", "w") as fo:
            fo.write("Index,sltbr_name,occupancy(%),frame/total,Distance(nm)±std.err\n")
            for i in range(len(sltbr_names)):
                name, occ, dist = sltbr_names[i], sltbr_occs[i], sltbr_dists[i]
                frame_num = np.sum(sltbr_issltbr[i])
                fo.write(f"{i},{name},{occ:.2%},{frame_num}/{t_num},{dist[0]}±{dist[1]}\n")

        ## output saltbridge number vs time to xvg file
        sltbr_number = np.sum(sltbr_issltbr, axis=0)
        xvg = XVG("SaltBridge_Number.xvg", new_file=True)
        xvg.title = "Salt Bridges Number"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Number"
        xvg.data_heads = ["Number"]
        xvg.data_columns = [time_array.tolist()] + [sltbr_number.tolist()]
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP SaltBridge module"
        xvg.save("SaltBridge_Number.xvg")
        cmd = f"""dit xvg_show -f SaltBridge_Number.xvg -ns -x "Time(ns)" -xs 0.001 -o SaltBridge_Number.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        # output distances to xvg
        xvg = XVG("SaltBridge_Distances.xvg", new_file=True)
        xvg.title = "Salt Bridges"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Distance(nm)"
        xvg.legends = sltbr_names
        xvg.data_heads = sltbr_names
        xvg.data_columns = [time_array.tolist()] + sltbr_datas
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP SaltBridge module"
        xvg.save("SaltBridge_Distances.xvg")
        cmd = f"""dit xvg_show -f SaltBridge_Distances.xvg -ns -x "Time(ns)" -xs 0.001 -o SaltBridge_Distance.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        # output issltbr to xpm
        xpm = XPM(f"SaltBridge_Existence_Map.xpm", new_file=True)
        xpm.height = len(sltbr_issltbr)
        xpm.width = t_num
        xpm.value_matrix = sltbr_issltbr[:]  ## deep copy
        xpm.title = "SaltBridge Existence Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "SaltBridge Index"
        xpm.type = "Discrete"
        xpm.xaxis = time_array.tolist()
        xpm.yaxis = [i for i in range(len(sltbr_names))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, sltbr_names)])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "SaltBridge"]
        xpm.colors = ["#FFFFFF", "#38A7D0"]
        xpm.chars = ["_", "o"]
        xpm.color_num = 2
        xpm.char_per_pixel = 1
        if len(xpm.dot_matrix) == 0 or len(xpm.datalines) == 0:
            for h in range(xpm.height):
                xpm.dot_matrix.append(["" for _ in range(xpm.width)])
                xpm.datalines.append("")
        for h in range(xpm.height):
            dot_line: str = ""
            for w in range(xpm.width):
                dot = xpm.chars[int(xpm.value_matrix[h][w])]
                xpm.dot_matrix[h][w] = dot
                dot_line += dot
            xpm.datalines[h] = dot_line
        xpm.save("SaltBridge_Existence_Map.xpm")
        cmd = f"""dit xpm_show -f SaltBridge_Existence_Map.xpm -ns -x "Time(ns)" -xs 0.001 -o SaltBridge_Existence_Map.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## remove black row or columns
        row2delete, column2delete = [], []
        for x in range(x_num):
            if np.sum(sltbr_occ_map[x, :]) == 0:
                row2delete.append(x)
        for y in range(y_num):
            if np.sum(sltbr_occ_map[:, y]) == 0:
                column2delete.append(y)
        sltbr_occ_map = np.delete(sltbr_occ_map, row2delete, axis=0)
        sltbr_occ_map = np.delete(sltbr_occ_map, column2delete, axis=1)
        sltbr_dist_ave_map = np.delete(sltbr_dist_ave_map, row2delete, axis=0)
        sltbr_dist_ave_map = np.delete(sltbr_dist_ave_map, column2delete, axis=1)
        sltbr_dist_std_map = np.delete(sltbr_dist_std_map, row2delete, axis=0)
        sltbr_dist_std_map = np.delete(sltbr_dist_std_map, column2delete, axis=1)
        positive_residues = np.delete(positive_residues, row2delete, axis=0)
        negative_residues = np.delete(negative_residues, column2delete, axis=0)
        x_num, y_num = sltbr_occ_map.shape

        # draw matrix for occupancy and distances
        sltbr_occ_map *= 100
        occ_mask = np.ones((x_num, y_num), dtype=float) - (sltbr_occ_map > 0)
        plt.clf()
        plt.figure(figsize=(15, 10))
        sns.heatmap(
            data=sltbr_occ_map,
            mask=occ_mask,
            vmin=0,
            vmax=100,
            center=50,
            cmap="Blues",
            cbar=True,
            cbar_kws={"label": "Occupancy(%)"},
            linewidth=0.1,
            annot=True,
            fmt=".2f",
        )
        plt.ylabel("Positive AA")
        plt.xlabel("Negative AA")
        plt.title("Occupancy(%)")
        plt.gca().invert_yaxis()
        plt.yticks([i + 0.5 for i in range(x_num)], positive_residues, rotation=0)
        plt.xticks([i + 0.5 for i in range(y_num)], negative_residues, rotation=90)
        plt.tight_layout()
        plt.savefig(f"""SaltBridges_Occupancy_Matrix.{self.conf["fig"]}""", dpi=300)
        plt.close()

        dist_mask = np.ones((x_num, y_num), dtype=float) - (sltbr_dist_ave_map > 0)
        plt.clf()
        plt.figure(figsize=(15, 10))
        sns.heatmap(
            data=sltbr_dist_ave_map,
            mask=dist_mask,
            vmin=0,
            vmax=dist_cutoff,
            cmap="Blues",
            cbar=True,
            cbar_kws={"label": "Distance Average(nm)"},
            linewidth=0.1,
            annot=True,
            fmt=".3f",
        )
        plt.ylabel("Positive AA")
        plt.xlabel("Negative AA")
        plt.title("Distance Average(nm)")
        plt.gca().invert_yaxis()
        plt.yticks([i + 0.5 for i in range(x_num)], positive_residues, rotation=0)
        plt.xticks([i + 0.5 for i in range(y_num)], negative_residues, rotation=90)
        plt.tight_layout()
        plt.savefig(f"""SaltBridges_Distance_Average_Matrix.{self.conf["fig"]}""", dpi=300)
        plt.close()

        plt.clf()
        plt.figure(figsize=(15, 10))
        sns.heatmap(
            data=sltbr_dist_std_map,
            mask=dist_mask,
            vmin=0,
            cmap="Blues",
            cbar=True,
            cbar_kws={"label": "Distance Std.Err(nm)"},
            linewidth=0.1,
            annot=True,
            fmt=".4f",
        )
        plt.ylabel("Positive AA")
        plt.xlabel("Negative AA")
        plt.title("Distance Std.Err(nm)")
        plt.gca().invert_yaxis()
        plt.yticks([i + 0.5 for i in range(x_num)], positive_residues, rotation=0)
        plt.xticks([i + 0.5 for i in range(y_num)], negative_residues, rotation=90)
        plt.tight_layout()
        plt.savefig(f"""SaltBridges_Distance_StdErr_Matrix.{self.conf["fig"]}""", dpi=300)
        plt.close()

        if not calc_lifetime:
            return 

        lifetime_curves = []
        for id, issltbr in enumerate(sltbr_issltbr):
            data = [set([int(d)]) if d == True else set() for d in issltbr]
            intermittent_data = correct_intermittency(data, intermittency)
            taus, curve, _ = autocorrelation(intermittent_data, tau_max, window_step)
            lifetime_curves.append(curve)
        if len(lifetime_curves) == 0:
            return 
        lifetime_curves = np.array(lifetime_curves)
        lt_time = np.array(taus) * (time_array[1]- time_array[0])
        xvg = XVG(f"SaltBridge_lifetime.xvg", new_file=True)
        xvg.title = "SaltBridge lifetime C(tau)"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "C(tau)"
        xvg.legends = sltbr_names
        xvg.data_heads = sltbr_names
        xvg.data_columns = [lt_time.tolist()] + lifetime_curves.tolist()
        xvg.row_num = len(lt_time)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP SaltBridge module"
        xvg.save(f"SaltBridge_lifetime.xvg")
        cmd = f"""dit xvg_show -f SaltBridge_lifetime.xvg -ns -o SaltBridge_lifetime.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        lifetimes = []
        for hbl_data in lifetime_curves:
            simps_values = simpson(hbl_data, lt_time)
            lifetimes.append(simps_values) # ps, lifetime by integration

        # output info to csv
        with open("SaltBridge_info.csv", "w") as fo:
            fo.write("Index,sltbr_name,occupancy(%),frame/total,lifetime(ps),Distance(nm)±std.err\n")
            for i in range(len(sltbr_names)):
                name, occ, dist = sltbr_names[i], sltbr_occs[i], sltbr_dists[i]
                lifetime = lifetimes[i]
                frame_num = np.sum(sltbr_issltbr[i])
                fo.write(f"{i},{name},{occ:.2%},{frame_num}/{t_num},{lifetime:.2f},{dist[0]:.4f}±{dist[1]:.4f}\n")
