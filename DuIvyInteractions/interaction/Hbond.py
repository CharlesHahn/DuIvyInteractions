"""
This module is part of DuIvyProcedures.procedures, designed for dealing Hbonds. 
Written by 杜艾维.
"""

import os
import sys
import numpy as np
from scipy.integrate import simpson
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBbondAna
from matplotlib import pyplot as plt
from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM


base = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
if base not in sys.path:
    sys.path.insert(0, base)

from utils import log
from framework.confParser import Config


class Hbond(log):
    def __init__(self, config: Config) -> None:
        self.conf = config
        self.load_mplstyle()

    @log.module_decorator
    def __call__(self) -> None:

        tpr = self.conf["tpr"]
        xtc = self.conf["xtc"]
        donor_group = self.conf["donor_group"]
        acceptor_group = self.conf["acceptor_group"]
        update_selection = self.conf["update_selection"]
        d_h_cutoff = self.conf["d_h_cutoff"]*10.0 # convert nm to A
        d_a_cutoff = self.conf["d_a_cutoff"]*10.0 # convert nm to A
        d_h_a_angle_cutoff = self.conf["d_h_a_angle_cutoff"]
        only_calc_number = self.conf["only_calc_number"]
        top2show = self.conf["top2show"]
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
        hb = HBbondAna(u, d_h_cutoff=d_h_cutoff, d_a_cutoff=d_a_cutoff, d_h_a_angle_cutoff=d_h_a_angle_cutoff, update_selections=update_selection)
        hb.hydrogens_sel = hb.guess_hydrogens(donor_group)
        hb.acceptors_sel = hb.guess_acceptors(acceptor_group)
        hb.run(start=fstart, stop=fend, step=fstep)
        
        ## count hbond number by frame
        xvg = XVG("hbnum.xvg", new_file=True)
        xvg.title = "hbond number"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Count"
        xvg.data_columns = [hb.times, hb.count_by_time()]
        xvg.data_heads = ["Count"]
        xvg.row_num = len(hb.times)
        xvg.column_num = 2
        xvg.save("hbnum.xvg")
        cmd = f"""dit xvg_show -f hbnum.xvg -ns -o hbnum.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)
        if only_calc_number:
            return 

        ## count hbond number by existance
        d = u.atoms[hb.results.hbonds[:, 1].astype(np.intp)]
        h = u.atoms[hb.results.hbonds[:, 2].astype(np.intp)]
        a = u.atoms[hb.results.hbonds[:, 3].astype(np.intp)]
        dha_list = np.array([d.ids, h.ids, a.ids]).T
        hbond_ids, ids_counts = np.unique(dha_list, axis=0, return_counts=True)
        unique_hbonds = np.concatenate((hbond_ids, ids_counts[:, None]), axis=1)
        unique_hbonds = unique_hbonds[unique_hbonds[:, 3].argsort()[::-1]]

        hbond_indexs = [f"{d}-{h}-{a}" for d, h, a in unique_hbonds[:, :3]]
        hbond_occupancy = [count/len(hb.times) for count in unique_hbonds[:, 3]]
        hbond_matrix = np.zeros((unique_hbonds.shape[0], len(hb.times)))
        hbond_angles = [[] for _ in range(unique_hbonds.shape[0])]
        hbond_distances = [[] for _ in range(unique_hbonds.shape[0])]
        hbond_names = ["" for _ in range(unique_hbonds.shape[0])]
        for hbidx, (donor, hydrogen, acceptor) in enumerate(unique_hbonds[:, :3]):
            da = u.atoms[int(donor)]
            d_name = f"{da.resname}{da.resnum}{da.name}({int(donor)+1})"
            ha = u.atoms[int(hydrogen)]
            h_name = f"{ha.resname}{ha.resnum}{ha.name}({int(hydrogen)+1})"
            aa = u.atoms[int(acceptor)]
            a_name = f"{aa.resname}{aa.resnum}{aa.name}({int(acceptor)+1})"
            hbond_names[hbidx] = f"{d_name}@{h_name}...{a_name}"

        ## NOTE: to avoid trajectory time start not from 0
        hbtimes = hb.times - u.trajectory[0].time
        frame_list = (hbtimes/u.trajectory.dt).astype(int).tolist()
        frame_index = {value:idx for idx, value in enumerate(frame_list)}
        hbond_index_dict = {value:idx for idx, value in enumerate(hbond_indexs)}
        for frame, donor, hydrogen, acceptor, distance, angle in hb.results.hbonds:
            hbidx = hbond_index_dict[f"{int(donor)}-{int(hydrogen)}-{int(acceptor)}"]
            fidx = frame_index[int(frame)]
            hbond_matrix[hbidx, fidx] += 1
            hbond_angles[hbidx].append(angle)
            hbond_distances[hbidx].append(distance*0.1) # convert to nm
        
        hbond_angle_ave = [np.average(angles) for angles in hbond_angles]
        hbond_distance_ave = [np.average(distances) for distances in hbond_distances]
        hbond_angle_std = [np.std(angles, ddof=1) for angles in hbond_angles]
        hbond_distance_std = [np.std(distances, ddof=1) for distances in hbond_distances]

        ## write all hbonds data
        with open("HBONDS_data.csv", "w") as fo:
            fo.write("id,donor@hydrogen...acceptor,occupancy(%),Present/Frames")
            fo.write(",Distance Ave(nm),Distance Std.err(nm),")
            fo.write("Angle Ave(deg),Angle Std.err(deg)\n")
            all_frame_num = len(hb.times)
            for id, count in enumerate(unique_hbonds[:, 3]):
                fo.write(f"{id},{hbond_names[id]},{hbond_occupancy[id]*100.0:.2f},")
                fo.write(f"{count:>d}/{all_frame_num:<d}")
                fo.write(f",{hbond_distance_ave[id]:>6.4f}")
                fo.write(f",{hbond_distance_std[id]:<6.4f}")
                fo.write(f",{hbond_angle_ave[id]:>6.2f}")
                fo.write(f",{hbond_angle_std[id]:<6.2f}\n")
        
        ## output hbond_matrix to xpm file
        xpm = XPM(f"hbmap.xpm", new_file=True)
        xpm.height = len(hbond_names)
        xpm.width = all_frame_num
        xpm.value_matrix = hbond_matrix.tolist()
        xpm.title = "Hbond Existence Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "Hbond Index"
        xpm.type = "Discrete"
        xpm.xaxis = hb.times.tolist()
        xpm.yaxis = [i for i in range(len(hbond_names))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, hbond_names)])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "Hbond"]
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
        xpm.save("hbmap.xpm")
        cmd = f"""dit xpm_show -f hbmap.xpm -ns -x "Time(ns)" -xs 0.001 -o hbmap.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        xpm = XPM(f"hbmap_top{top2show}.xpm", new_file=True)
        xpm.height = len(hbond_names[:top2show])
        xpm.width = all_frame_num
        xpm.value_matrix = hbond_matrix.tolist()[:top2show]
        xpm.title = "Hbond Existence Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "Hbond Index"
        xpm.type = "Discrete"
        xpm.xaxis = hb.times.tolist()
        xpm.yaxis = [i for i in range(len(hbond_names[:top2show]))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, hbond_names[:top2show])])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "Hbond"]
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
        xpm.save(f"hbmap_top{top2show}.xpm")
        cmd = f"""dit xpm_show -f hbmap_top{top2show}.xpm -ns -x "Time(ns)" -xs 0.001 -o hbmap_top{top2show}.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)
        
        ## calculate hbond distribution vs amino acids
        donor_distribution = {}
        acceptor_distribution = {}
        for frame, donor, hydrogen, acceptor, distance, angle in hb.results.hbonds:
            da = u.atoms[int(donor)]
            d_name = f"{da.resname}{da.resnum}"
            aa = u.atoms[int(acceptor)]
            a_name = f"{aa.resname}{aa.resnum}"
            if d_name not in donor_distribution:
                donor_distribution[d_name] = 1
            else:
                donor_distribution[d_name] += 1
            if a_name not in acceptor_distribution:
                acceptor_distribution[a_name] = 1
            else:
                acceptor_distribution[a_name] += 1

        donor_distribution = sorted(donor_distribution.items(), key=lambda x:x[1], reverse=True)
        acceptor_distribution = sorted(acceptor_distribution.items(), key=lambda x:x[1], reverse=True)
        with open("hbond_distribution.dat", 'w') as fo:
            fo.write("Donor, Count\n")
            for d, c in donor_distribution:
                fo.write(f"{d:<10s}{c:>5d}\n")
            fo.write("\n")
            fo.write("Acceptor, Count\n")
            for a, c in acceptor_distribution:
                fo.write(f"{a:<10s}{c:>5d}\n")
            fo.write("\n")
        
        ## draw the distribution of hbond number
        plt.clf()
        plt.bar([d[0] for d in donor_distribution[:top2show]], [d[1] for d in donor_distribution[:top2show]])
        plt.xticks(rotation=90)
        plt.xlabel("Amino Acid")
        plt.ylabel("Hbond Count")
        plt.title(f"""Top {top2show} hbond distribution by donor""")
        plt.tight_layout()
        plt.savefig(f"hbond_donor_distribution.{self.conf['fig']}")
        plt.close()
        plt.clf()
        plt.bar([a[0] for a in acceptor_distribution[:top2show]], [a[1] for a in acceptor_distribution[:top2show]])
        plt.xticks(rotation=90)
        plt.xlabel("Amino Acid")
        plt.ylabel("Hbond Count")
        plt.title(f"""Top {top2show} hbond distribution by acceptor""")
        plt.tight_layout()
        plt.savefig(f"hbond_acceptor_distribution.{self.conf['fig']}")
        plt.close()

        if not calc_lifetime: 
            return 

        hbond_lifetime_curves = []
        for id, (_, hidx, aidx) in enumerate(unique_hbonds[:top2show, :3]):
            hbl = HBbondAna(u, d_h_cutoff=d_h_cutoff, d_a_cutoff=d_a_cutoff, d_h_a_angle_cutoff=d_h_a_angle_cutoff, update_selections=False)
            hbl.hydrogens_sel = f"index {hidx}"
            hbl.acceptors_sel = f"index {aidx}"
            hbl.run(start=fstart, stop=fend, step=fstep)

            taus, hblife_curve, = hbl.lifetime(tau_max=tau_max, window_step=window_step, intermittency=intermittency)
            hbond_lifetime_curves.append(hblife_curve)
        if len(hbond_lifetime_curves) == 0:
            return 
        hbond_lifetime_curves = np.array(hbond_lifetime_curves)
        hbl_time = taus * (hb.times[1]- hb.times[0])
        xvg = XVG(f"hbond_top{top2show}_lifetime.xvg", new_file=True)
        xvg.title = "Hbond lifetime C(tau)"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "C(tau)"
        xvg.legends = hbond_names[:top2show]
        xvg.data_heads = hbond_names[:top2show]
        xvg.data_columns = [hbl_time.tolist()] + hbond_lifetime_curves.tolist()
        xvg.row_num = len(hbl_time)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP Hbond module"
        xvg.save(f"hbond_top{top2show}_lifetime.xvg")
        cmd = f"""dit xvg_show -f hbond_top{top2show}_lifetime.xvg -ns -o hbond_top{top2show}_lifetime.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        hbond_lifetimes = []
        for hbl_data in hbond_lifetime_curves:
            simps_values = simpson(hbl_data, hbl_time)
            hbond_lifetimes.append(simps_values) # ps, lifetime by integration

        ## write top2show hbonds data
        with open(f"HBONDS_data_top{top2show}.csv", "w") as fo:
            fo.write("id,donor@hydrogen...acceptor,occupancy(%),Present/Frames")
            fo.write(",Lifetime(ps)")
            fo.write(",Distance Ave(nm),Distance Std.err(nm),")
            fo.write("Angle Ave(deg),Angle Std.err(deg)\n")
            all_frame_num = len(hb.times)
            for id, count in enumerate(unique_hbonds[:top2show, 3]):
                fo.write(f"{id},{hbond_names[id]},{hbond_occupancy[id]*100.0:.2f},")
                fo.write(f"{count:>d}/{all_frame_num:<d}")
                fo.write(f",{hbond_lifetimes[id]:>6.4f}")
                fo.write(f",{hbond_distance_ave[id]:>6.4f}")
                fo.write(f",{hbond_distance_std[id]:<6.4f}")
                fo.write(f",{hbond_angle_ave[id]:>6.2f}")
                fo.write(f",{hbond_angle_std[id]:<6.2f}\n")