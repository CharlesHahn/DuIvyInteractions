"""
This module is part of DuIvyProcedures.procedures, designed for dealing PiCation. 
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


class PiCation(log):
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

    def calc_offset(self, coors, center_vec, dist):
        ring_vec0 = coors[2] - coors[0]
        ring_vec1 = coors[2] - coors[4]
        ring_normal = np.cross(ring_vec0, ring_vec1)
        deg_c = self.calc_degree(ring_normal, center_vec)
        offset = math.sin(deg_c / 180.0 * math.pi) * dist
        return offset

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
        group4PiRing = self.conf["group4PiRing"]
        group4Cation = self.conf["group4Cation"]
        byIndex = self.conf["byIndex"]
        only_aromatic_rings = self.conf["only_aromatic_rings"]
        other_ring_max_atom_num = self.conf["other_ring_max_atom_num"]
        planarity_cutoff = self.conf["planarity_cutoff"]
        Pi_rings_Index = self.conf["Pi_rings_Index"]
        Pi_rings_Index = np.array(Pi_rings_Index) - 1  # mda index start from 0

        NH3_atomnames = self.conf["NH3_atomnames"]
        COO_atomnames = self.conf["COO_atomnames"]
        Backbone_atomnames = self.conf[
            "Backbone_atomnames"
        ]  #  ["H", "N", "CA", "C", "O"]
        ignore_chain_end = self.conf["ignore_chain_end"]
        cation_Index = self.conf["cation_Index"]
        cation_Index = np.array(cation_Index) - 1  # mda atom index (start from 0)
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
                all_atoms.write("PiRings_byIndex.pdb")
            else:
                self.error("No Pi rings found in your yaml input !!!")
                return
            Pi_rings_names = []
            for id, ring_Index in enumerate(Pi_rings_Index):
                res = u.atoms[ring_Index].residues
                name = f"{res.resnames[0]}{res.resnums[0]}"
                if name in Pi_rings_names:
                    name += f"_{id}"
                Pi_rings_names.append(name)

            if len(cation_Index) > 0:
                all_atoms = u.atoms[[i for i in chain(*cation_Index)]]
                all_atoms.write("Cations_byIndex.pdb")
            else:
                self.error("No Cations found in your yaml input !!!")
                return
            cation_names = []
            for id, index in enumerate(cation_Index):
                res = u.atoms[index].residues
                name = f"{res.resnames[0]}{res.resnums[0]}"
                if name in cation_names:
                    name += f"_{id}"
                cation_names.append(name)

        else:  # by rdkit to find rings and cations
            elements = [mda.topology.guessers.guess_atom_element(n) for n in u.atoms.names]
            u.add_TopologyAttr("elements", elements)
            not_water = " and ".join([f"not resname {w}" for w in ["HOH", "WAT", "SOL"]])

            ## Pi rings
            Pi_atoms = u.select_atoms(f"{not_water} and {group4PiRing}")
            Pi_atoms.write("group4PiRing_4rdkit2findrings.pdb")
            Pi_aromatic_rings, Pi_other_rings = self.find_rings("group4PiRing_4rdkit2findrings.pdb")
            Pi_aromatic_rings = [sorted(list(r)) for r in Pi_aromatic_rings]
            Pi_other_rings = [
                sorted(list(r))
                for r in Pi_other_rings
                if len(r) <= other_ring_max_atom_num and len(r) >= 5
            ]
            Pi_aromatic_rings = [[Pi_atoms.ids[id] for id in r] for r in Pi_aromatic_rings]
            Pi_other_rings = [[Pi_atoms.ids[id] for id in r] for r in Pi_other_rings]
            planar_rings = []
            for other_ring in Pi_other_rings:    # check planarity
                coors = u.atoms[other_ring].positions * 0.1  # A to nm
                if self.calc_ring_planarity(coors, planarity_cutoff):
                    planar_rings.append(other_ring)
            Pi_other_rings = planar_rings
            self.info(
                f"Found {len(Pi_aromatic_rings)} aromatic rings and {len(Pi_other_rings)} other rings in group1 {group4PiRing}"
            )
            if len(Pi_aromatic_rings) > 0:
                u.atoms[[i for i in chain(*Pi_aromatic_rings)]].write(
                    "PiRings_Aromatic.pdb"
                )
            if len(Pi_other_rings) > 0:
                u.atoms[[i for i in chain(*Pi_other_rings)]].write(
                    "PiRings_Other.pdb"
                )
            if only_aromatic_rings == True:
                Pi_rings_Index = Pi_aromatic_rings
            else:
                Pi_rings_Index = Pi_aromatic_rings + Pi_other_rings
            if len(Pi_rings_Index) == 0:
                self.error(f"No Pi rings found in your structure of group4PiRing {group4PiRing} !!!")
                return
            Pi_rings_names = []
            for id, ring_Index in enumerate(Pi_rings_Index):
                res = u.atoms[ring_Index].residues
                name = f"{res.resnames[0]}{res.resnums[0]}"
                if name in Pi_rings_names:
                    name += f"_{id}"
                Pi_rings_names.append(name)

            ## group4Cation
            cation_uAAs, cation_residues = [], []
            cation_uAAc, cation_residues_c = [], []
            for res in u.select_atoms(group4Cation).residues:
                res_set = set(res.atoms.names)
                iswhole = set(Backbone_atomnames).issubset(res_set)
                if iswhole == True:  ## to obtain charge center of sidechain
                    key = " and ".join([f"not name {a}" for a in Backbone_atomnames])
                    side = res.atoms.select_atoms(key)
                    if np.sum(side.charges) > 0.42:
                        cation_uAAs.append(side)
                        cation_residues.append(f"{res.resname}{res.resnum}")
                else:  ## not whole, get C-ter, N-ter, and sidechains
                    if set(NH3_atomnames).issubset(res_set):
                        key = " or ".join([f"name {a}" for a in NH3_atomnames])
                        key += " or name CA"  # add CA, which always charged
                        NH3_atoms = res.atoms.select_atoms(key)
                        NH3_name = f"{res.resname}{res.resnum}_NH3"
                        if np.sum(NH3_atoms.charges) > 0.42:
                            cation_uAAc.append(NH3_atoms)
                            cation_residues_c.append(NH3_name)
                        else:
                            self.warn(
                                f"The NH3 group of {NH3_name} contained charges less than 0.42, quite wired? DIP ignored it."
                            )
                            print(NH3_atoms.charges)
                    NH3_COO_Backbone = set(
                        NH3_atomnames + COO_atomnames + Backbone_atomnames
                    )
                    key = " and ".join([f"not name {a}" for a in NH3_COO_Backbone])
                    side = res.atoms.select_atoms(key)
                    if np.sum(side.charges) > 0.42:
                        cation_uAAc.append(side)
                        cation_residues_c.append(f"{res.resname}{res.resnum}_side")
            self.info(
                f"Found {len(cation_uAAs)} AA with POSITIVE sidechain in system"
            )
            self.info(
                f"And found {len(cation_uAAc)} AA (broken backbone) with POSITIVE charge in system"
            )
            if not ignore_chain_end:
                cation_uAAs += cation_uAAc
                cation_residues += cation_residues_c

            cation_Index = [side.ids for side in cation_uAAs]
            cation_names = cation_residues
            if len(cation_Index) != 0:
                u.atoms[[i for i in chain(*cation_Index)]].write("Cations_AtomGroups.pdb")
            else:
                self.error("Empty Indexs, No charged groups found in your structure.")
                return

        ## output names of rings and cations, and indexs for checking and reuse this module
        with open("PiCation_Names_Indexs.txt", "w") as fo:
            fo.write(f"PiRings_Names, Indexs\n")
            for name, ring_Index in zip(Pi_rings_names, Pi_rings_Index):
                fo.write(f"{name}, {[r+1 for r in ring_Index]}\n")
            fo.write(f"Cations_Names, Indexs\n")
            for name, ring_Index in zip(cation_names, cation_Index):
                fo.write(f"{name}, {[r+1 for r in ring_Index]}\n")

        ## get the distance of ring centers
        time_array, dist_dataframe = [], []
        for ts in u.trajectory[fstart:fend:fstep]:
            Pi_rings_centers, cation_centers = [], []
            for ring_Index in Pi_rings_Index:
                Pi_rings_centers.append(u.atoms[ring_Index].center_of_mass())
            for ring_Index in cation_Index:
                cation_centers.append(u.atoms[ring_Index].center_of_charge())
            Pi_rings_centers = np.array(Pi_rings_centers)
            cation_centers = np.array(cation_centers)
            dist_matrix = mda_dist.distance_array(
                Pi_rings_centers, cation_centers, box=ts.dimensions
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
        new_XY, new_dist_dataframe, isPiCation_dataframe = [], [], []
        for x in range(x_num):
            for y in range(y_num):
                data = dist_OK_dataframe[:, x, y]
                if np.any(data):
                    new_XY.append((x, y))
                    new_dist_dataframe.append(dist_dataframe[:, x, y])
                    isPiCation_dataframe.append(data)
        new_dist_dataframe = np.array(new_dist_dataframe)
        isPiCation_dataframe = np.array(isPiCation_dataframe)

        ## do the filter by offset
        new_offset_dataframe = np.zeros(new_dist_dataframe.shape)
        for i, (x, y) in enumerate(new_XY):
            x_atom_indexs, y_atom_indexs = Pi_rings_Index[x], cation_Index[y]
            for t in range(t_num):
                ts = u.trajectory[fstart:fend:fstep][t]
                x_atom_coors = u.atoms[x_atom_indexs].positions
                x_center = u.atoms[x_atom_indexs].center_of_mass()
                y_center = u.atoms[y_atom_indexs].center_of_charge()
                offset = self.calc_offset(
                    x_atom_coors,
                    x_center - y_center,
                    new_dist_dataframe[i, t] * 10,
                )
                offset *= 0.1  # A to nm
                new_offset_dataframe[i, t] = offset
                if offset > ring_center_offset:  ## offset criteria
                    isPiCation_dataframe[i, t] = False
            self.info(
                f"Calculating {i}/{len(new_XY)} PiStacking for {Pi_rings_names[x]}-{cation_names[y]}......"
            )

        PiCation_names = [f"{Pi_rings_names[x]}-{cation_names[y]}" for x, y in new_XY]
        PiCation_occupancy, pication2delete = [], []
        for i, (x, y) in enumerate(new_XY):
            occupancy = np.sum(isPiCation_dataframe[i, :]) / t_num
            PiCation_occupancy.append(occupancy)
            if occupancy == 0:
                pication2delete.append(i)

        ## to delete blank PiCation
        new_XY = np.delete(new_XY, pication2delete, axis=0)
        PiCation_names = np.delete(PiCation_names, pication2delete, axis=0)
        new_dist_dataframe = np.delete(new_dist_dataframe, pication2delete, axis=0)
        new_offset_dataframe = np.delete(new_offset_dataframe, pication2delete, axis=0)
        isPiCation_dataframe = np.delete(isPiCation_dataframe, pication2delete, axis=0)
        PiCation_occupancy = np.delete(PiCation_occupancy, pication2delete, axis=0)
        ## to avoild NO PiCation
        if len(new_XY) == 0:
            self.error("No PiCation found in your structure !!!")
            return

        ## output the information to csv file
        with open("PiCation_Info.csv", "w") as fo:
            fo.write(
                f"id,Name,Occupancy,Frames/Total,Distance(nm),Offset(nm)\n"
            )
            for i, (x, y) in enumerate(new_XY):
                dist = np.mean(
                    new_dist_dataframe[i, :][isPiCation_dataframe[i, :]], axis=0
                )
                offset = np.mean(
                    new_offset_dataframe[i, :][isPiCation_dataframe[i, :]], axis=0
                )
                fo.write(
                    f"{i},{PiCation_names[i]},{PiCation_occupancy[i]:.2%},{np.sum(isPiCation_dataframe[i, :])}/{t_num},{dist:.6f},{offset:.6f}\n"
                )

        ## output pication number vs time to xvg file
        pication_number = np.sum(isPiCation_dataframe, axis=0)
        xvg = XVG("PiCation_Number.xvg", new_file=True)
        xvg.title = "PiCation Number"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Number"
        xvg.data_heads = ["Number"]
        xvg.data_columns = [time_array.tolist()] + [pication_number.tolist()]
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiCation module"
        xvg.save("PiCation_Number.xvg")
        cmd = f"""dit xvg_show -f PiCation_Number.xvg -ns -x "Time(ns)" -xs 0.001 -o PiCation_Number.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## output pication distance and offset vs time to xvg file
        xvg = XVG("PiCation_Distances.xvg", new_file=True)
        xvg.title = "PiCation Distances"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Distances(nm)"
        xvg.legends = PiCation_names
        xvg.data_heads = PiCation_names
        xvg.data_columns = [time_array.tolist()] + new_dist_dataframe.tolist()
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiCation module"
        xvg.save("PiCation_Distances.xvg")
        cmd = f"""dit xvg_show -f PiCation_Distances.xvg -ns -x "Time(ns)" -xs 0.001 -o PiCation_Distances.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)
        xvg = XVG("PiCation_Offsets.xvg", new_file=True)
        xvg.title = "PiCation Offsets"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "Offsets(nm)"
        xvg.legends = PiCation_names
        xvg.data_heads = PiCation_names
        xvg.data_columns = [time_array.tolist()] + new_offset_dataframe.tolist()
        xvg.row_num = len(time_array)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiCation module"
        xvg.save("PiCation_Offsets.xvg")
        cmd = f"""dit xvg_show -f PiCation_Offsets.xvg -ns -x "Time(ns)" -xs 0.001 -o PiCation_Offsets.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        ## output isPiCation to xpm file
        xpm = XPM(f"PiCation_Existence_Map.xpm", new_file=True)
        xpm.height = len(PiCation_names)
        xpm.width = t_num
        xpm.value_matrix = isPiCation_dataframe.tolist()
        xpm.title = "PiCation Existence Map"
        xpm.xlabel = "Time(ps)"
        xpm.ylabel = "PiCation Index"
        xpm.type = "Discrete"
        xpm.xaxis = time_array.tolist()
        xpm.yaxis = [i for i in range(len(PiCation_names))]
        xpm.legend = " ".join([f"{i}:{n}" for i, n in zip(xpm.yaxis, PiCation_names)])
        xpm.yaxis.reverse()  ## xpm store data and yaxis from high to low
        xpm.value_matrix.reverse()
        xpm.notes = ["None", "PiCation"]
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
        xpm.save("PiCation_Existence_Map.xpm")
        cmd = f"""dit xpm_show -f PiCation_Existence_Map.xpm -ns -x "Time(ns)" -xs 0.001 -o PiCation_Existence_Map.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        if not calc_lifetime:
            return 

        lifetime_curves = []
        for id, isPiCation in enumerate(isPiCation_dataframe):
            data = [set([int(d)]) if d == True else set() for d in isPiCation]
            intermittent_data = correct_intermittency(data, intermittency)
            taus, curve, _ = autocorrelation(intermittent_data, tau_max, window_step)
            lifetime_curves.append(curve)
        if len(lifetime_curves) == 0:
            return 
        lifetime_curves = np.array(lifetime_curves)
        lt_time = np.array(taus) * (time_array[1]- time_array[0])
        xvg = XVG(f"PiCation_lifetime.xvg", new_file=True)
        xvg.title = "PiCation lifetime C(tau)"
        xvg.xlabel = "Time(ps)"
        xvg.ylabel = "C(tau)"
        xvg.legends = PiCation_names
        xvg.data_heads = PiCation_names
        xvg.data_columns = [lt_time.tolist()] + lifetime_curves.tolist()
        xvg.row_num = len(lt_time)
        xvg.column_num = len(xvg.data_columns)
        xvg.comments = "## generated by DIP PiCation module"
        xvg.save(f"PiCation_lifetime.xvg")
        cmd = f"""dit xvg_show -f PiCation_lifetime.xvg -ns -o PiCation_lifetime.{self.conf["fig"]}"""
        status, output, error = self.run_terminal(cmd)

        lifetimes = []
        for hbl_data in lifetime_curves:
            simps_values = simpson(hbl_data, lt_time)
            lifetimes.append(simps_values) # ps, lifetime by integration

        ## output the information to csv file
        with open("PiCation_Info.csv", "w") as fo:
            fo.write(
                f"id,Name,Occupancy,Frames/Total,lifetime(ps),Distance(nm),Offset(nm)\n"
            )
            for i, (x, y) in enumerate(new_XY):
                dist = np.mean(
                    new_dist_dataframe[i, :][isPiCation_dataframe[i, :]], axis=0
                )
                offset = np.mean(
                    new_offset_dataframe[i, :][isPiCation_dataframe[i, :]], axis=0
                )
                fo.write(
                    f"{i},{PiCation_names[i]},{PiCation_occupancy[i]:.2%},{np.sum(isPiCation_dataframe[i, :])}/{t_num},{lifetimes[i]:.2f},{dist:.6f},{offset:.6f}\n"
                )

