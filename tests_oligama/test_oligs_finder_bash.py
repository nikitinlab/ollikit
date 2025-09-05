from ollikit import AddOligsNew
import os
import tempfile
import pytest
import numpy as np

def test_add_oligs_new_basic():
    """Тест базового функционала AddOligsNew"""
    input_data = {
        "target_seqs": ["ACTGCTAGAGATTTTCCACAT"],
        "target_names": ["LucS8-B1"],
        "sequence_types": ["DNA"],
        "add_name": "D99",
        "min_w0x0": 0.9,
        "max_w0x0": 0.995,
        "bad_thr": 0.02,
        "rounds": 30,
        "olig_conc": 1e-6,
        "celsius": 25,
        "num_complex": 2,
        "bad_energy": -0.1,
        "metric": "fraction",
        "aff_predictor": "Nupack",
        "hairpin_predictor": "Nupack"
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        finder = AddOligsNew(input_data, tmpdir)
        result = finder.run(timeout_seconds=60, min_required=1)
        
        # Проверяем структуру результата
        assert "results" in result
        assert "log" in result
        assert isinstance(result["results"], list)
        assert isinstance(result["log"], str)
        
        # Если нашли кандидатов
        if result["results"]:
            candidate = result["results"][0]
            assert "name" in candidate
            assert "seq" in candidate
            assert "energy" in candidate
            assert "affinity_to_target" in candidate
            assert "pattern" in candidate
            
            # Проверяем значения
            assert candidate["name"] == "D99"
            assert len(candidate["seq"]) == len(input_data["target_seqs"][0])
            assert isinstance(candidate["energy"], float)
            assert isinstance(candidate["affinity_to_target"], float)
            assert isinstance(candidate["pattern"], str)

def test_add_oligs_new_validation():
    """Тест валидации входных параметров"""
    # Тест с некорректными порогами
    with pytest.raises(Exception):
        input_data = {
            "target_seqs": ["ACTGCTAGAGATTTTCCACAT"],
            "target_names": ["LucS8-B1"],
            "sequence_types": ["DNA"],
            "min_w0x0": 0.99,  # min > max
            "max_w0x0": 0.9,
            "metric": "fraction"
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            AddOligsNew(input_data, tmpdir)

    # Тест с множественными таргетами
    with pytest.raises(Exception):
        input_data = {
            "target_seqs": ["ACTGCTAGAGATTTTCCACAT", "TGTGGAAAATCTCTAGCAGTT"],
            "target_names": ["seq1", "seq2"],
            "sequence_types": ["DNA", "DNA"],
            "metric": "fraction"
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            AddOligsNew(input_data, tmpdir)

def test_add_oligs_new_defaults():
    """Тест значений по умолчанию"""
    input_data = {
        "target_seqs": ["ACTGCTAGAGATTTTCCACAT"],
        "target_names": ["LucS8-B1"],
        "sequence_types": ["DNA"],
        "metric": "fraction"
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        finder = AddOligsNew(input_data, tmpdir)
        
        # Проверяем значения по умолчанию
        assert finder.add_name == "D99"
        assert finder.min_w0x0 == 0.9
        assert finder.max_w0x0 == 0.995
        assert finder.bad_thr == 0.02
        assert finder.rounds == 30
        assert finder.olig_conc == 1e-6
        assert finder.num_complex == 2
        assert finder.bad_energy == -0.1

def test_add_oligs_new_pattern():
    """Тест генерации паттерна"""
    input_data = {
        "target_seqs": ["ACTGCTAGAGATTTTCCACAT"],
        "target_names": ["LucS8-B1"],
        "sequence_types": ["DNA"],
        "metric": "fraction"
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        finder = AddOligsNew(input_data, tmpdir)
        
        # Тест _build_pattern
        pattern = finder._build_pattern("ACTG", "ACTG")
        assert pattern == "ACTG"
        
        pattern = finder._build_pattern("ACTG", "ACTT")
        assert pattern == "ACTx"
        
        pattern = finder._build_pattern("ACT", "ACTG")
        assert pattern == "ACTx"

def test_add_oligs_new_complement():
    """Тест комплементарных последовательностей"""
    input_data = {
        "target_seqs": ["ACTGCTAGAGATTTTCCACAT"],
        "target_names": ["LucS8-B1"],
        "sequence_types": ["DNA"],
        "metric": "fraction"
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        finder = AddOligsNew(input_data, tmpdir)
        
        # Тест _reverse_complement
        assert finder._reverse_complement("ACTG") == "CAGT"
        assert finder._reverse_complement("AAAA") == "TTTT"
        assert finder._reverse_complement("GCTA") == "TAGC"

