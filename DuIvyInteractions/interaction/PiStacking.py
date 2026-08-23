"""
This module is part of DuIvyProcedures.procedures, designed for dealing PiStacking. 
Written by 杜艾维.
"""

import os
import sys
import math
import numpy as np
import MDAnalysis as mda
from itertools import chain
from scipy.integrate import simpson
from MDAnalysis.analysis import distances as mda_dist
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
from DuIvyTools.DuIvyTools.FileParser.xvgParser import XVG
from DuIvyTools.DuIvyTools.FileParser.xpmParser import XPM

base = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
if base not in sys.path:
    sys.path.insert(0, base)

from utils import log
from framework.confParser import Config


class PiStacking(log):
    def __init__(self, config: Config) -> None:
        self.conf = config

    def find_rings(self, pdbfile):
        from rdkit import Chem

        mol = Chem.MolFromPDBFile(pdbfile, removeHs=False)
        if mol == None:
            self.error("Error in reading pdb by rdkit in calc PiStacking")
            return [], []
        rings = mol.GetRingInfo()
        if not rings:
            self.error("No ring in your structure")
            return [], []
        aromatic_rings, other_rings = [], []
        bond_rings = rings.BondRings()
        atom_rings = rings.AtomRings()
        for atom_ring, bond_ring in zip(atom_rings, bond_rings):
            bond_prop = [mol.GetBondWithIdx(idx).GetIsAromatic() for idx in bond_ring]
            if np.all(bond_prop):
                aromatic_rings.append(atom_ring)
            else:
                other_rings.append(atom_ring)
        return aromatic_rings, other_rings

    def calc_degree(self, vec1, vec2):
        dotProduct = np.dot(vec1, vec2)
        vec1_length = np.linalg.norm(vec1)
        vec2_length = np.linalg.norm(vec2)
        cos_degree = dotProduct / (vec1_length * vec2_length)
        if cos_degree > 1:
            cos_degree = 1.0
        elif cos_degree < -1:
            cos_degree = -1.0
        try:
            degree = math.acos(cos_degree) * 180 / math.pi
        except:
            self.critical(f"crash in calc degree: cos_deg= {cos_degree}")
        degree = math.acos(cos_degree) * 180 / math.pi
        if degree > 90:
            degree = 180 - degree
        return degree

    def calc_pistacking(self, coors1, coors2, center_vec, dist):
        ring1_vec0 = coors1[2] - coors1[0]
        ring1_vec1 = coors1[2] - coors1[4]
        ring1_normal = np.cross(ring1_vec0, ring1_vec1)
        ring2_vec0 = coors2[2] - coors2[0]
        ring2_vec1 = coors2[2] - coors2[4]
        ring2_normal = np.cross(ring2_vec0, ring2_vec1)
        deg = self.calc_degree(ring1_normal, ring2_normal)
        deg_1c = self.calc_degree(ring1_normal, center_vec)
        offset_1 = math.sin(deg_1c / 180.0 * math.pi) * dist
        deg_2c = self.calc_degree(ring2_normal, center_vec)
        offset_2 = math.sin(deg_2c / 180.0 * math.pi) * dist
        offset = np.min([offset_1, offset_2])
        return deg, offset

    def calc_ring_planarity(self, coors, planarity_cutoff):
        fake_ind = [i for i in range(len(coors))] * 2
        for i in range(len(coors)):
            for j in range(i + 1, len(coors)):
                i_vec0 = coors[fake_ind[i+1]] - coors[fake_ind[i]]
                i_vec1 = coors[fake_ind[i+1]] - coors[fake_ind[i+2]]
                i_normal = np.cross(i_vec0, i_vec1)
                j_vec0 = coors[fake_ind[j+1]] - coors[fake_ind[j]]
                j_vec1 = coors[fake_ind[j+1]] - coors[fake_ind[j+2]]
                j_normal = np.cross(j_vec0, j_vec1)
                deg = self.calc_degree(i_normal, j_normal)
                if deg > planarity_cutoff:
                    return False
        return True

    @log.module_decorator
    def __call__(self) -> None:

        tpr = self.conf["tpr"]
        xtc = self.conf["xtc"]
        dist_max_cutoff = self.conf["distance_max_cutoff"]
        dist_min_cutoff = self.conf["distance_min_cutoff"]
        ring_center_offset = self.conf["ring_center_offset"]
        angle_T_min, angle_T_max = self.conf["angle4T_stacking"]
        angle_P_min, angle_P_max = self.conf["angle4P_stacking"]
        group1 = self.conf["group1"]
        group2 = self.conf["group2"]
        byIndex = self.conf["byIndex"]
        only_aromatic_rings = self.conf["only_aromatic_rings"]
        other_ring_max_atom_num = self.conf["other_ring_max_atom_num"]
        planarity_cutoff = self.conf["planarity_cutoff"]
        Pi_rings_Index = self.conf["Pi_rings_Index"]
        Pi_rings_Index = np.array(Pi_rings_Index) - 1  # mda index start from 0
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
            if len(Pi_rings_Index) > 0:
                all_atoms = u.atoms[[i for i in chain(*Pi_rings_Index)]]
                all_atoms.write("PiStack_byIndex.pdb")
            if len(Pi_rings_Index) == 0:
                self.error("No Pi rings found in your structure !!!")
                return
            g1_Pi_rings_Index = Pi_rings_Index
            g2_Pi_rings_Index = Pi_rings_Index
        else:  # by rdkit to find rings
            elements = [mda.topology.guessers.guess_atom_element(n) for n in u.atoms.names]
            u.add_TopologyAttr("elements", elements)
            not_water = " and ".join([f"not resname {w}" for w in ["HOH", "WAT", "SOL"]])

            ## group1
            g1_atoms = u.select_atoms(f"{not_water} and {group1}")
            g1_atoms.write("group1_allatoms4rdkit2findrings.pdb")
            g1_aromatic_rings, g1_other_rings = self.find_rings("group1_allatoms4rdkit2findrings.pdb")
            g1_aromatic_rings = [sorted(list(r)) for r in g1_aromatic_rings]
            g1_other_rings = [
                sorted(list(r))
                for r in g1_other_rings
                if len(r) <= other_ring_max_atom_num and len(r) >= 5
            ]
            g1_aromatic_rings = [[g1_atoms.ids[id] for id in r] for r in g1_aromatic_rings]
            g1_other_rings = [[g1_atoms.ids[id] for id in r] for r in g1_other_rings]
            planar_rings = []
            for other_ring in g1_other_rings:    # check planarity
                coors = u.atoms[other_ring].positions * 0.1  # A to nm
                if self.calc_ring_planarity(coors, planarity_cutoff):
                    planar_rings.append(other_ring)
            g1_other_rings = planar_rings
            self.info(
                f"Found {len(g1_aromatic_rings)} aromatic rings and {len(g1_other_rings)} other rings in group1 {group1}"
            )
            if len(g1_aromatic_rings) > 0:
                u.atoms[[i for i in chain(*g1_aromatic_rings)]].write(
                    "g1_PiStack_byRings_Aromatic.pdb"
                )
            if len(g1_other_rings) > 0:
                u.atoms[[i for i in chain(*g1_other_rings)]].write(
                    "g1_PiStack_byRings_Other.pdb"
                )
            if only_aromatic_rings == True:
                g1_Pi_rings_Index = g1_aromatic_rings
            else:
                g1_Pi_rings_Index = g1_aromatic_rings + g1_other_rings

            ## group2
            g2_atoms = u.select_atoms(f"{not_water} and {group2}")
            g2_atoms.write("group2_allatoms4rdkit2findrings.pdb")
            g2_aromatic_rings, g2_other_rings = self.find_rings("group2_allatoms4rdkit2findrings.pdb")
            g2_aromatic_rings = [sorted(list(r)) for r in g2_aromatic_rings]
            g2_other_rings = [
                sorted(list(r))
                for r in g2_other_rings
                if len(r) <= other_ring_max_atom_num and len(r) >= 5
            ]
            g2_aromatic_rings = [[g2_atoms.ids[id] for id in r] for r in g2_aromatic_rings]
            g2_other_rings = [[g2_atoms.ids[id] for id in r] for r in g2_other_rings]
            planar_rings = []
            for other_ring in g2_other_rings:    # check planarity
                coors = u.atoms[other_ring].positions * 0.1  # A to nm
                if self.calc_ring_planarity(coors, planarity_cutoff):
                    planar_rings.append(other_ring)
            g2_other_rings = planar_rings
            self.info(
                f"Found {len(g2_aromatic_rings)} aromatic rings and {len(g2_other_rings)} other rings in group2 {group2}"
            )
            if len(g2_aromatic_rings) > 0:
                u.atoms[[i for i in chain(*g2_aromatic_rings)]].write(
                    "g2_PiStack_byRings_Aromatic.pdb"
                )
            if len(g2_other_rings) > 0:
                u.atoms[[i for i in chain(*g2_other_rings)]].write(
                    "g2_PiStack_byRings_Other.pdb"
                )
            if only_aromatic_rings == True:
                g2_Pi_rings_Index = g2_aromatic_rings
            else:
                g2_Pi_rings_Index = g2_aromatic_rings + g2_other_rings

            if len(g1_Pi_rings_Index) == 0:
                self.error(f"No Pi rings found in your structure of group1 {group1} !!!")
                return
            if len(g2_Pi_rings_Index) == 0:
                self.error(f"No Pi rings found in your structure of group2 {group2} !!!")
                return

        ## get the ben ring names
        g1_ben_names, g2_ben_names = [], []
        for id, ring_Index in enumerate(g1_Pi_rings_Index):
            res = u.atoms[ring_Index].residues
            name = f"{res.resnames[0]}{res.resnums[0]}"
            if name in g1_ben_names:
                name += f"_{id}"
            g1_ben_names.append(name)
        for id, ring_Index in enumerate(g2_Pi_rings_Index):
            res = u.atoms[ring_Index].residues
            name = f"{res.resnames[0]}{res.resnums[0]}"
            if name in g2_ben_names:
                name += f"_{id}"
            g2_ben_names.append(name)
        ## output rings names and indexs for checking and reuse this module
        with open("PiStacking_Names_Indexs.txt", "w") as fo:
            fo.write(f"Group1 PiStacking_Names, Indexs\n")
            for name, ring_Index in zip(g1_ben_names, g1_Pi_rings_Index):
                fo.write(f"{name}, {[r+1 for r in ring_Index]}\n")
            fo.write(f"Group2 PiStacking_Names, Indexs\n")
            for name, ring_Index in zip(g2_ben_names, g2_Pi_rings_Index):
                fo.write(f"{name}, {[r+1 for r in ring_Index]}\n")

        ## get the distance of ring centers
        time_array, dist_dataframe = [], []
        for ts in u.trajectory[fstart:fend:fstep]:
            g1_ben_centers, g2_ben_centers = [], []
            for ring_Index in g1_Pi_rings_Index:
                g1_ben_centers.append(u.atoms[ring_Index].center_of_mass())
            for ring_Index in g2_Pi_rings_Index:
                g2_ben_centers.append(u.atoms[ring_Index].center_of_mass())
            g1_ben_centers = np.array(g1_ben_centers)
            g2_ben_centers = np.array(g2_ben_centers)
            dist_matrix = mda_dist.distance_array(
                g1_ben_centers, g2_ben_centers, box=ts.dimensions
            )
            dist_matrix *= 0.1  # A to nm
            time_array.append(ts.time)
            dist_dataframe.append(dist_matrix)
        time_array = np.array(time_array)
        dist_dataframe = np.array(dist_dataframe)

        # do the filter by center distance cutoff
        t_num, x_num, y_num = dist_dataframe.shape
        dist_OK_dataframe = np.logical_and(
            (dist_dataframe <= dist_max_cutoff), (dist_dataframe >= dist_min_cutoff)
        )  ## distance criteria
        new_XY, new_dist_dataframe, isPistack_dataframe = [], [], []
        for x in range(x_num):
            for y in range(y_num):
                data = dist_OK_dataframe[:, x, y]
                if np.any(data):
                    new_XY.append((x, y))
                    new_dist_dataframe.append(dist_dataframe[:, x, y])
                    isPistack_dataframe.append(data)
        new_dist_dataframe = np.array(new_dist_dataframe)
        isPistack_dataframe = np.array(isPistack_dataframe)

        ## to delete neighbor ring pairs
        delete_neighbor_rings = []
        for i, (x, y) in enumerate(new_XY):
            if len(set(g1_Pi_rings_Index[x]) & set(g2_Pi_rings_Index[y])) > 0:
                print(g1_Pi_rings_Index[x], g2_Pi_rings_Index[y])
                delete_neighbor_rings.append(i)
        new_XY = np.delete(new_XY, delete_neighbor_rings, axis=0)
        new_dist_dataframe = np.delete(
            new_dist_dataframe, delete_neighbor_rings, axis=0
        )
        isPistack_dataframe = np.delete(
            isPistack_dataframe, delete_neighbor_rings, axis=0
        )
        ## to delete repeat ring pairs
        all_names, delete_repeat_rings = [], []
        for i, (x, y) in enumerate(new_XY):
            i_name = f"{g1_ben_names[x]}-{g2_ben_names[y]}"
            repeat_name = f"{g2_ben_names[y]}-{g1_ben_names[x]}"
            if repeat_name in all_names:
                delete_repeat_rings.append(i)
                print(f"{i_name} repeated to {repeat_name}, delete it !")
            else:
                all_names.append(i_name)
        new_XY = np.delete(new_XY, delete_repeat_rings, axis=0)
        new_dist_dataframe = np.delete(
            new_dist_dataframe, delete_repeat_rings, axis=0
        )
        isPistack_dataframe = np.delete(
            isPistack_dataframe, delete_repeat_rings, axis=0
        )

        ## do the filter by angles
        new_angle_dataframe = np.zeros(new_dist_dataframe.shape)
        new_offset_dataframe = np.zeros(new_dist_dataframe.shape)
        for i, (x, y) in enumerate(new_XY):
            x_atom_indexs, y_atom_indexs = g1_Pi_rings_Index[x], g2_Pi_rings_Index[y]
            for t in range(t_num):
                ts = u.trajectory[fstart:fend:fstep][t]
                x_atom_coors = u.atoms[x_atom_indexs].positions
                y_atom_coors = u.atoms[y_atom_indexs].positions
                x_center = u.atoms[x_atom_indexs].center_of_mass()
                y_center = u.atoms[y_atom_indexs].center_of_mass()
                angle, offset = self.calc_pistacking(
                    x_atom_coors,
                    y_atom_coors,
                    x_center - y_center,
                    new_dist_dataframe[i, t] * 10,
                )
                offset *= 0.1  # A to nm
                new_angle_dataframe[i, t] = angle
                new_offset_dataframe[i, t] = offset
                if offset > ring_center_offset:  ## offset criteria
                    isPistack_dataframe[i, t] = False
                if angle > angle_P_max and angle < angle_T_min:  ## angle criteria
                    isPistack_dataframe[i, t] = False
            self.info(
                f"Calculating {i}/{len(new_XY)} PiStacking for {g1_ben_names[x]}-{g2_ben_names[y]}......"
            )

        PiStack_names = [f"{g1_ben_names[x]}-{g2_ben_names[y]}" for x, y in new_XY]
        pistack_type_dataframe = np.zeros(new_dist_dataframe.shape)
        pistack_occupancy, pistack2delete = [], []
        for i, (x, y) in enumerate(new_XY):
            isPistack_data = isPistack_dataframe[i, :]
            occupancy = np.sum(isPistack_data) / t_num
            pistack_occupancy.append(occupancy)
            if occupancy == 0:
                pistack2delete.append(i)
            else:
                is_indexs = np.where(isPistack_data)[0]
                for index in is_indexs:
                    angle = new_angle_dataframe[i, index]
                    if angle >= angle_P_min and angle <= angle_P_max:
                        pistack_type_dataframe[i, index] = 1  # P-Stacking
                    elif angle >= angle_T_min and angle <= angle_T_max:
                        pistack_type_dataframe[i, index] = 2  # T-Stacking

        ## to delete blank PiStacking
        new_XY = np.delete(new_XY, pistack2delete, axis=0)
        PiStack_names = np.delete(PiStack_names, pistack2delete, axis=0)
        new_dist_dataframe = np.delete(new_dist_dataframe, pistack2delete, axis=0)
        new_angle_dataframe = np.delete(new_angle_dataframe, pistack2delete, axis=0)
        new_offset_dataframe = np.delete(new_offset_dataframe, pistack2delete, axis=0)
        isPistack_dataframe = np.delete(isPistack_dataframe, pistack2delete, axis=0)
        pistack_occupancy = np.delete(pistack_occupancy, pistack2delete, axis=0)
        pistack_type_dataframe = np.delete(
            pistack_type_dataframe, pistack2delete, axis=0
        )
        ## to avoild NO pistacking
        if len(new_XY) == 0:
            self.error("No PiStacking found in your structure !!!")
            return

        ## output the information to csv file
        with open("PiStacking_Info.csv", "w") as fo:
            fo.write(
                f"id,Name,Occupancy,Distance(nm),Offset(nm),P-Stacking_Occupancy,T-Stacking_Occupancy,P-Angle(deg),T-Angle(deg),Total_Frames,P-Frames,T-Frames\n"
            )
            for i, (x, y) in enumerate(new_XY):
                dist = np.mean(
                    new_dist_dataframe[i, :][isPistack_dataframe[i, :]], axis=0
                )
                offset = np.mean(
                    new_offset_dataframe[i, :][isPistack_dataframe[i, :]], axis=0
                )
                P_indexs = pistack_type_dataframe[i, :] == 1
                T_indexs = pistack_type_dataframe[i, :] == 2
                P_occ = np.sum(P_indexs) / t_num
                T_occ = np.sum(T_indexs) / t_num
                P_angle = np.mean(new_angle_dataframe[i, :][P_indexs], axis=0)
                T_angle = np.mean(new_angle_dataframe[i, :][T_indexs], axis=0)
                fo.write(
                    f"{i},{PiStack_names[i]},{pistack_occupancy[i]:.2%},{dist:.6f},{offset:.6f},{P_occ:.2%},{T_occ:.2%},{P_angle:.2f},{T_angle:.2f},{t_num},{np.sum(P_indexs)},{np.sum(T_indexs)}\n"
                )

        ## output pistacking number vs time to xvg file
        pistacking_number = np.sum(isPistack_dataframe, axis=0)
        P_pistacking_number = np.sum(pistack_type_dataframe==1, axis=0)
        T_pistacking_number = np.sum(pistack_type_dataframe==2, axis=0)
        xvg = XVG("PiStacking_Number.xvg", new_file=True)
        xvg.title = "PiStacking Number"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Number"
        xvg.legends = ["Total", "P-Stacking", "T-Stacking"]
        xvg.data_heads = ["Total", "P-Stacking", "T-Stacking"]
        xvg.data_columns = [time_array.tolist()] + [pistacking_number.tolist()]
        xvg.data_columns += [P_pistacking_number.tolist()]
        xvg.data_columns += [T_pistacking_number.tolist()]
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiStacking module"
        xvg.save("PiStacking_Number.xvg")
        cmd = f"""dit xvg_show -f PiStacking_Number.xvg -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Number.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## output pistacking distances, angles, and offsets to xvg file
        xvg = XVG("PiStacking_Distances.xvg", new_file=True)
        xvg.title = "PiStacking Distances"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Distances(nm)"
        xvg.legends = PiStack_names
        xvg.data_heads = PiStack_names
        xvg.data_columns = [time_array.tolist()] + new_dist_dataframe.tolist()
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiStacking module"
        xvg.save("PiStacking_Distances.xvg")
        cmd = f"""dit xvg_show -f PiStacking_Distances.xvg -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Distances.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)
        xvg = XVG("PiStacking_Angles.xvg", new_file=True)
        xvg.title = "PiStacking Angles"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Angle(deg)"
        xvg.legends = PiStack_names
        xvg.data_heads = PiStack_names
        xvg.data_columns = [time_array.tolist()] + new_angle_dataframe.tolist()
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiStacking module"
        xvg.save("PiStacking_Angles.xvg")
        cmd = f"""dit xvg_show -f PiStacking_Angles.xvg -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Angles.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)
        xvg = XVG("PiStacking_Offsets.xvg", new_file=True)
        xvg.title = "PiStacking Offsets"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Offsets(nm)"
        xvg.legends = PiStack_names
        xvg.data_heads = PiStack_names
        xvg.data_columns = [time_array.tolist()] + new_offset_dataframe.tolist()
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiStacking module"
        xvg.save("PiStacking_Offsets.xvg")
        cmd = f"""dit xvg_show -f PiStacking_Offsets.xvg -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Offsets.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## output isPistack to xpm file
        xpm = XPM(f"PiStacking_Existence_Map.xpm", new_file=True)
        xpm.height = len(PiStack_names)
        xpm.width = t_num
        xpm.value_matrix = isPistack_dataframe.tolist()
        xpm.title = "PiStacking Existence Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "PiStacking Index"
        xpm.type = "Discrete"
        xpm.xaxis = time_array.tolist()
        xpm.yaxis = [i for i in range(len(PiStack_names))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, PiStack_names)])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "PiStacking"]
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
        xpm.save("PiStacking_Existence_Map.xpm")
        cmd = f"""dit xpm_show -f PiStacking_Existence_Map.xpm -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Existence_Map.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## output PiStack_type to xpm file
        xpm = XPM(f"PiStacking_Type_Map.xpm", new_file=True)
        xpm.height = len(PiStack_names)
        xpm.width = t_num
        xpm.value_matrix = pistack_type_dataframe.tolist()
        xpm.title = "PiStacking Type Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "PiStacking Index"
        xpm.type = "Discrete"
        xpm.xaxis = time_array.tolist()
        xpm.yaxis = [i for i in range(len(PiStack_names))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, PiStack_names)])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "P-Stacking", "T-Stacking"]
        xpm.colors = ["#FFFFFF", "#38A7D0", "#F67088"]
        xpm.chars = ["_", "o", "x"]
        xpm.color_num = 3
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
        xpm.save("PiStacking_Type_Map.xpm")
        cmd = f"""dit xpm_show -f PiStacking_Type_Map.xpm -ns -x "Time(ns)" -xs 0.001 -o PiStacking_Type_Map.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        if not calc_lifetime:
            return 

        lifetime_curves = []
        for id, isPistack_data in enumerate(isPistack_dataframe):
            data = [set([int(d)]) if d == True else set() for d in isPistack_data]
            intermittent_data = correct_intermittency(data, intermittency)
            taus, curve, _ = autocorrelation(intermittent_data, tau_max, window_step)
            lifetime_curves.append(curve)
        if len(lifetime_curves) == 0:
            return 
        lifetime_curves = np.array(lifetime_curves)
        lt_time = np.array(taus) * (time_array[1]- time_array[0])
        xvg = XVG(f"PiStacking_lifetime.xvg", new_file=True)
        xvg.title = "PiStacking lifetime C(tau)"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "C(tau)"
        xvg.legends = PiStack_names
        xvg.data_heads = PiStack_names
        xvg.data_columns = [lt_time.tolist()] + lifetime_curves.tolist()
        xvg.row_num = len(lt_time)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiStacking module"
        xvg.save(f"PiStacking_lifetime.xvg")
        cmd = f"""dit xvg_show -f PiStacking_lifetime.xvg -ns -o PiStacking_lifetime.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        lifetimes = []
        for hbl_data in lifetime_curves:
            simps_values = simpson(hbl_data, lt_time)
            lifetimes.append(simps_values) # ps, lifetime by integration

        ## output the information to csv file
        with open("PiStacking_Info.csv", "w") as fo:
            fo.write(
                f"id,Name,Occupancy,Distance(nm),Offset(nm),lifetime(ps),P-Stacking_Occupancy,T-Stacking_Occupancy,P-Angle(deg),T-Angle(deg),Total_Frames,P-Frames,T-Frames\n"
            )
            for i, (x, y) in enumerate(new_XY):
                dist = np.mean(
                    new_dist_dataframe[i, :][isPistack_dataframe[i, :]], axis=0
                )
                offset = np.mean(
                    new_offset_dataframe[i, :][isPistack_dataframe[i, :]], axis=0
                )
                P_indexs = pistack_type_dataframe[i, :] == 1
                T_indexs = pistack_type_dataframe[i, :] == 2
                P_occ = np.sum(P_indexs) / t_num
                T_occ = np.sum(T_indexs) / t_num
                P_angle = np.mean(new_angle_dataframe[i, :][P_indexs], axis=0)
                T_angle = np.mean(new_angle_dataframe[i, :][T_indexs], axis=0)
                fo.write(
                    f"{i},{PiStack_names[i]},{pistack_occupancy[i]:.2%},{dist:.6f},{offset:.6f},{lifetimes[i]:.2f},{P_occ:.2%},{T_occ:.2%},{P_angle:.2f},{T_angle:.2f},{t_num},{np.sum(P_indexs)},{np.sum(T_indexs)}\n"
                )

