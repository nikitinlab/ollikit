import os
import random
import tempfile
from pathlib import Path

import pytest
from ollikit import OligsFinderBash


def test_oligs_finder_bash_basic():
    """
    Минимальный интеграционный тест с nupack.
    """
    random.seed(42)

    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        # Минимальный InputList: таргет + два «прочих» олига
        (root / "InputList").write_text(
            "LucS8-B1\tATGCGTACGATC\n"
            "S1\tTTTTTTTTTT\n"
            "S2\tCCCCCCCCCC\n",
            encoding="utf-8",
        )

        runner = OligsFinderBash(
            target_name="LucS8-B1",
            add_name="D99",
            project_root=root,
            persist_files=False,    
            material="dna",
            temperature_c=25,
            rounds=12,               # чтобы тест был быстрым
            min_w0x0=0.10,
            max_w0x0=0.95,
            bad_thr=0.05,
            bad_energy=-0.5,
        )

        out = runner.run(timeout_seconds=20, min_required=1)

        assert isinstance(out, dict)
        assert "results" in out and isinstance(out["results"], list)
        assert len(out["results"]) >= 1

        first = out["results"][0]
        # проверки структуры результата
        assert "name" in first and "seq" in first
        assert "energy" in first and isinstance(first["energy"], float)
        assert "affinity_to_target" in first and isinstance(first["affinity_to_target"], float)
        assert "pattern" in first

        # при persist_files=False должен вернуться  InputList
        assert out.get("updated_input_list_text")
        assert "D99" in out["updated_input_list_text"]
        assert first["seq"] in out["updated_input_list_text"]


def test_oligs_finder_bash_persist_files():
    """
    Вариант с persist_files=True: тут реально дописываем InputList и лог.
    """
    random.seed(123)

    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        input_list = root / "InputList"
        addtmp = root / "AddTMP"
        aux = root / "AUX"

        input_list.write_text(
            "LucS8-B1\tATGCGTACGATC\n"
            "S1\tTTTTTTTTTT\n"
            "S2\tCCCCCCCCCC\n",
            encoding="utf-8",
        )

        runner = OligsFinderBash(
            target_name="LucS8-B1",
            add_name="D99",
            project_root=root,
            persist_files=True,      # теперь пишем на диск
            material="dna",
            temperature_c=25,
            rounds=10,
            min_w0x0=0.10,
            max_w0x0=0.95,
            bad_thr=0.05,
            bad_energy=-0.5,
        )

        out = runner.run(timeout_seconds=20, min_required=1)

        # Есть результаты
        assert out["results"]

        # Проверяем, что строка реально дописалась в InputList
        updated_text = input_list.read_text(encoding="utf-8")
        assert "D99" in updated_text
        assert out["results"][0]["seq"] in updated_text

        # Лог и AUX должны существовать и быть не пустыми
        assert addtmp.exists()
        assert aux.exists()
        assert addtmp.read_text(encoding="utf-8").strip()
