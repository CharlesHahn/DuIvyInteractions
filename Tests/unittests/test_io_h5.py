# -*- coding: utf-8 -*-
"""HDF5 存储和加载的单元测试。

测试往返一致性：保存后加载，验证数据无损。
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path

from DuIvyInteractions.core.datas import Interaction, Group, AtomData
from DuIvyInteractions.io.h5 import save_interactions, load_interactions


# 临时目录路径
TEMP_DIR = Path(__file__).parent.parent.parent / "test_temp"


@pytest.fixture
def sample_atom():
    """创建示例 AtomData。"""
    return AtomData(
        atom_global_idx=100,
        atom_idx_in_residue=5,
        atom_name="CG",
        atom_type="ca",
        atom_element="C",
        atom_charge=-0.1,
        atom_mass=12.011
    )


@pytest.fixture
def sample_group(sample_atom):
    """创建示例 Group。"""
    return Group(
        group_id=1,
        group_type="aromatic_ring",
        molecule="D927",
        residue_name="D927",
        residue_id=200,
        atoms=[sample_atom],
        metadata={"source": "amber", "ring_size": 6}
    )


@pytest.fixture
def sample_interaction(sample_group):
    """创建示例 Interaction。"""
    n_pairs = 2
    n_frames = 10
    
    # 创建两个不同的 group
    group1 = sample_group
    group2 = Group(
        group_id=2,
        group_type="H_donor",
        molecule="RBD",
        residue_name="ARG",
        residue_id=73,
        atoms=[
            AtomData(
                atom_global_idx=200,
                atom_idx_in_residue=10,
                atom_name="NE",
                atom_type="N",
                atom_element="N",
                atom_charge=-0.5,
                atom_mass=14.007
            ),
            AtomData(
                atom_global_idx=201,
                atom_idx_in_residue=11,
                atom_name="HE",
                atom_type="H",
                atom_element="H",
                atom_charge=0.3,
                atom_mass=1.008
            )
        ],
        metadata={"source": "residue_name"}
    )
    
    return Interaction(
        interaction_type="hydrogen_bond",
        groups=[(group1, group2), (group2, group1)],
        existence=np.array([[True, False, True, False, True, False, True, False, True, False],
                            [False, True, False, True, False, True, False, True, False, True]]),
        metrics={
            "distance": np.random.rand(n_pairs, n_frames) * 4.0,
            "angle": np.random.rand(n_pairs, n_frames) * 180.0
        }
    )


@pytest.fixture
def empty_interaction():
    """创建空的 Interaction。"""
    return Interaction(
        interaction_type="pi_stacking",
        groups=[],
        existence=np.empty((0, 10), dtype=bool),
        metrics={"distance": np.empty((0, 10)), "angle": np.empty((0, 10))}
    )


@pytest.fixture(autouse=True)
def setup_temp_dir():
    """创建临时目录。"""
    TEMP_DIR.mkdir(exist_ok=True)
    yield
    # 测试后清理
    import shutil
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR, ignore_errors=True)


class TestH5Roundtrip:
    """测试往返一致性。"""

    def test_roundtrip_single_interaction(self, sample_interaction):
        """测试单个 Interaction 的往返。"""
        path = TEMP_DIR / "test.h5"
        
        # 保存
        save_interactions([sample_interaction], str(path))
        
        # 加载
        loaded = load_interactions(str(path))
        
        # 验证数量
        assert len(loaded) == 1
        
        # 验证基本信息
        loaded_interaction = loaded[0]
        assert loaded_interaction.interaction_type == sample_interaction.interaction_type
        assert loaded_interaction.n_pairs == sample_interaction.n_pairs
        assert loaded_interaction.n_frames == sample_interaction.n_frames
        
        # 验证 existence
        np.testing.assert_array_equal(loaded_interaction.existence, sample_interaction.existence)
        
        # 验证 metrics
        for key in sample_interaction.metrics:
            np.testing.assert_array_almost_equal(
                loaded_interaction.metrics[key],
                sample_interaction.metrics[key]
            )

    def test_roundtrip_groups_structure(self, sample_interaction):
        """测试 groups 结构的往返。"""
        path = TEMP_DIR / "test_groups.h5"
        
        save_interactions([sample_interaction], str(path))
        loaded = load_interactions(str(path))
        
        loaded_groups = loaded[0].groups
        original_groups = sample_interaction.groups
        
        # 验证 tuple 数量
        assert len(loaded_groups) == len(original_groups)
        
        # 验证每个 tuple
        for loaded_tuple, original_tuple in zip(loaded_groups, original_groups):
            assert len(loaded_tuple) == len(original_tuple)
            
            for loaded_group, original_group in zip(loaded_tuple, original_tuple):
                # 验证 group 基本信息
                assert loaded_group.group_id == original_group.group_id
                assert loaded_group.group_type == original_group.group_type
                assert loaded_group.molecule == original_group.molecule
                assert loaded_group.residue_name == original_group.residue_name
                assert loaded_group.residue_id == original_group.residue_id
                
                # 验证 atoms 数量
                assert len(loaded_group.atoms) == len(original_group.atoms)
                
                # 验证每个 atom
                for loaded_atom, original_atom in zip(loaded_group.atoms, original_group.atoms):
                    assert loaded_atom.atom_global_idx == original_atom.atom_global_idx
                    assert loaded_atom.atom_idx_in_residue == original_atom.atom_idx_in_residue
                    assert loaded_atom.atom_name == original_atom.atom_name
                    assert loaded_atom.atom_type == original_atom.atom_type
                    assert loaded_atom.atom_element == original_atom.atom_element
                    assert abs(loaded_atom.atom_charge - original_atom.atom_charge) < 1e-10
                    assert abs(loaded_atom.atom_mass - original_atom.atom_mass) < 1e-10

    def test_roundtrip_metadata(self, sample_interaction):
        """测试 metadata 的往返。"""
        path = TEMP_DIR / "test_metadata.h5"
        
        save_interactions([sample_interaction], str(path))
        loaded = load_interactions(str(path))
        
        loaded_groups = loaded[0].groups
        original_groups = sample_interaction.groups
        
        for loaded_tuple, original_tuple in zip(loaded_groups, original_groups):
            for loaded_group, original_group in zip(loaded_tuple, original_tuple):
                # 验证 metadata
                assert loaded_group.metadata == original_group.metadata

    def test_roundtrip_multiple_interactions(self, sample_interaction):
        """测试多个 Interaction 的往返。"""
        path = TEMP_DIR / "test_multiple.h5"
        
        # 创建第二个 interaction
        interaction2 = Interaction(
            interaction_type="pi_stacking",
            groups=sample_interaction.groups[:1],
            existence=sample_interaction.existence[:1],
            metrics={"distance": sample_interaction.metrics["distance"][:1]}
        )
        
        interactions = [sample_interaction, interaction2]
        
        save_interactions(interactions, str(path))
        loaded = load_interactions(str(path))
        
        # 验证数量
        assert len(loaded) == 2
        
        # 验证每个 interaction
        for loaded_int, original_int in zip(loaded, interactions):
            assert loaded_int.interaction_type == original_int.interaction_type
            assert loaded_int.n_pairs == original_int.n_pairs

    def test_roundtrip_empty_interaction(self, empty_interaction):
        """测试空 Interaction 的往返。"""
        path = TEMP_DIR / "test_empty.h5"
        
        save_interactions([empty_interaction], str(path))
        loaded = load_interactions(str(path))
        
        assert len(loaded) == 1
        assert loaded[0].interaction_type == "pi_stacking"
        assert loaded[0].n_pairs == 0
        assert loaded[0].n_frames == 10

    def test_roundtrip_with_none_metadata(self):
        """测试 metadata 包含 None 值的往返。"""
        atom = AtomData(
            atom_global_idx=0,
            atom_idx_in_residue=0,
            atom_name="C",
            atom_type="ca",
            atom_element="C",
            atom_charge=0.0,
            atom_mass=12.011
        )
        
        group = Group(
            group_id=1,
            group_type="aromatic_ring",
            molecule="D927",
            residue_name="D927",
            residue_id=1,
            atoms=[atom],
            metadata={"source": None, "ring_size": 6}
        )
        
        interaction = Interaction(
            interaction_type="pi_stacking",
            groups=[(group, group)],
            existence=np.array([[True, False]]),
            metrics={"distance": np.array([[3.5, 4.5]])}
        )
        
        path = TEMP_DIR / "test_none.h5"
        
        save_interactions([interaction], str(path))
        loaded = load_interactions(str(path))
        
        assert loaded[0].groups[0][0].metadata["source"] is None
        assert loaded[0].groups[0][0].metadata["ring_size"] == 6

    def test_roundtrip_with_compression(self, sample_interaction):
        """测试压缩模式的往返。"""
        path_compressed = TEMP_DIR / "test_compressed.h5"
        path_uncompressed = TEMP_DIR / "test_uncompressed.h5"
        
        save_interactions([sample_interaction], str(path_compressed), compress=True)
        save_interactions([sample_interaction], str(path_uncompressed), compress=False)
        
        # 验证两种模式都能正确加载
        loaded_compressed = load_interactions(str(path_compressed))
        loaded_uncompressed = load_interactions(str(path_uncompressed))
        
        assert loaded_compressed[0].interaction_type == loaded_uncompressed[0].interaction_type
        np.testing.assert_array_equal(
            loaded_compressed[0].existence,
            loaded_uncompressed[0].existence
        )


class TestH5EdgeCases:
    """测试边界情况。"""

    def test_large_interaction(self):
        """测试大型 Interaction。"""
        n_pairs = 100
        n_frames = 1000
        
        groups = []
        for i in range(n_pairs):
            atom = AtomData(
                atom_global_idx=i,
                atom_idx_in_residue=0,
                atom_name="C",
                atom_type="ca",
                atom_element="C",
                atom_charge=-0.1,
                atom_mass=12.011
            )
            group = Group(
                group_id=i,
                group_type="aromatic_ring",
                molecule="D927",
                residue_name="D927",
                residue_id=i,
                atoms=[atom],
                metadata={"index": i}
            )
            groups.append((group, group))
        
        interaction = Interaction(
            interaction_type="pi_stacking",
            groups=groups,
            existence=np.random.rand(n_pairs, n_frames) > 0.5,
            metrics={
                "distance": np.random.rand(n_pairs, n_frames) * 5.0,
                "angle": np.random.rand(n_pairs, n_frames) * 90.0
            }
        )
        
        path = TEMP_DIR / "test_large.h5"
        
        save_interactions([interaction], str(path))
        loaded = load_interactions(str(path))
        
        assert loaded[0].n_pairs == n_pairs
        assert loaded[0].n_frames == n_frames
        assert len(loaded[0].groups) == n_pairs

    def test_special_characters_in_strings(self):
        """测试字符串中的特殊字符。"""
        atom = AtomData(
            atom_global_idx=0,
            atom_idx_in_residue=0,
            atom_name="Cα",
            atom_type="ca",
            atom_element="C",
            atom_charge=0.0,
            atom_mass=12.011
        )
        
        group = Group(
            group_id=1,
            group_type="aromatic_ring",
            molecule="D927-Ligand",
            residue_name="D927",
            residue_id=1,
            atoms=[atom],
            metadata={"description": "π-π stacking between aromatic rings"}
        )
        
        interaction = Interaction(
            interaction_type="pi_stacking",
            groups=[(group, group)],
            existence=np.array([[True]]),
            metrics={"distance": np.array([[3.5]])}
        )
        
        path = TEMP_DIR / "test_special.h5"
        
        save_interactions([interaction], str(path))
        loaded = load_interactions(str(path))
        
        assert loaded[0].groups[0][0].molecule == "D927-Ligand"
        assert loaded[0].groups[0][0].metadata["description"] == "π-π stacking between aromatic rings"

    def test_tuple_length_varies(self):
        """测试不同长度的 tuple。"""
        atom = AtomData(
            atom_global_idx=0,
            atom_idx_in_residue=0,
            atom_name="C",
            atom_type="ca",
            atom_element="C",
            atom_charge=0.0,
            atom_mass=12.011
        )
        
        group1 = Group(
            group_id=1,
            group_type="H_donor",
            molecule="D927",
            residue_name="D927",
            residue_id=1,
            atoms=[atom],
            metadata={}
        )
        
        group2 = Group(
            group_id=2,
            group_type="H_acceptor",
            molecule="D927",
            residue_name="D927",
            residue_id=2,
            atoms=[atom],
            metadata={}
        )
        
        group3 = Group(
            group_id=3,
            group_type="water",
            molecule="SOL",
            residue_name="SOL",
            residue_id=3,
            atoms=[atom],
            metadata={}
        )
        
        # 2-tuple 和 3-tuple 混合
        interaction = Interaction(
            interaction_type="water_bridge",
            groups=[
                (group1, group2),      # 2-tuple
                (group1, group3, group2)  # 3-tuple
            ],
            existence=np.array([[True, False], [False, True]]),
            metrics={"distance": np.array([[3.5, 4.5], [5.5, 6.5]])}
        )
        
        path = TEMP_DIR / "test_tuple.h5"
        
        save_interactions([interaction], str(path))
        loaded = load_interactions(str(path))
        
        assert len(loaded[0].groups) == 2
        assert len(loaded[0].groups[0]) == 2
        assert len(loaded[0].groups[1]) == 3
