"""Run an existing small input with both references and the fixed right tail."""
import csv
import pathlib
import subprocess
import sys
import tempfile

root = pathlib.Path(__file__).resolve().parents[1]
binary = root / "build/scgsca_run"
input_dir = pathlib.Path(sys.argv[1]).resolve()
for flags in (["--ct-tail"], ["--ct-tail", "right"], ["--ct-tail", "left"], ["--ct-tail=left"]):
    result = subprocess.run([str(binary), *flags], capture_output=True, text=True, timeout=10)
    assert result.returncode != 0 and "ct-tail" in result.stderr, result
help_result = subprocess.run([str(binary), "--help"], capture_output=True, text=True, check=True)
assert "fixed right tail" in help_result.stdout
assert "Permutation.Pvalue.AbsDelta" in help_result.stdout
assert "--ct-tail" not in help_result.stdout

with tempfile.TemporaryDirectory(prefix="ct_fixed_tail_") as tmp:
    rows = {}
    for ref in ("A", "B"):
        out = pathlib.Path(tmp) / ref
        cmd = [str(binary), "--kegg-dir", str(input_dir), "--out-dir", str(out),
               "--files", "1", "--nperm", "5", "--cores", "1", "--seed", "416627",
               "--gene-start", "3", "--ref-letter", ref, "--unpaired", "--skip-dc",
               "--ct-pair", "Epithelial cell,Immune"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        files = list(out.rglob("*_delta_result.tsv"))
        assert len(files) == 1, result.stdout
        with files[0].open() as handle:
            rows[ref] = next(csv.DictReader(handle, delimiter="\t"))
        assert rows[ref]["CT.Tail"] == "right"
        assert float(rows[ref]["Nperm.Paired"]) == 5
        other = "B" if ref == "A" else "A"
        assert f"CT tail: right; Delta = D({other}) - D({ref})" in result.stdout
    assert abs(float(rows["A"]["Observed.Delta"])+float(rows["B"]["Observed.Delta"])) < 1e-5
    # The unchanged v5 RNG can produce different permutation samples after relabeling.
    for ref in rows:
        assert 1/6 - 1e-6 <= float(rows[ref]["Permutation.Pvalue.Delta"]) <= 1
        assert 1/6 - 1e-6 <= float(rows[ref]["Permutation.Pvalue.AbsDelta"]) <= 1
        print(ref, rows[ref])
print("CLI and end-to-end fixed-right-tail smoke tests passed.")
