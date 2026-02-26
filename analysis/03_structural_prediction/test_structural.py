"""test_structural.py -- Unit tests for Step 3: Structural Prediction.

Run: pytest analysis/03_structural_prediction/ -v
"""

import importlib
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add script directory to path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

mod01 = importlib.import_module("01_select_candidates")
mod02 = importlib.import_module("02_fetch_sequences")
mod03 = importlib.import_module("03_generate_lsf_jobs")
mod04 = importlib.import_module("04_parse_results")

PROJECT_ROOT = SCRIPT_DIR.parent.parent
CONSENSUS_PATH = PROJECT_ROOT / "results/01_cross_cohort_overlap/step01_consensus_proteins_pos.csv"
ENRICHMENT_DIR = PROJECT_ROOT / "results/02_pathway_enrichment"


# --- Test candidate selection ---

class TestSelectCandidates:

    @pytest.fixture
    def candidates(self):
        if not CONSENSUS_PATH.exists():
            pytest.skip("Step 1 output not found")
        return mod01.select_candidates(CONSENSUS_PATH, ENRICHMENT_DIR, n_candidates=20)

    def test_vegfa_included(self, candidates):
        assert "VEGFA" in candidates["target"].values, "VEGFA positive control missing"

    def test_positive_controls_present(self, candidates):
        controls = {"VEGFA", "NRP1", "NRP2", "SEMA3A"}
        present = set(candidates["target"].values)
        missing = controls - present
        # Controls might already be in consensus; check they appear somewhere
        assert len(missing) <= 1, f"Missing controls: {missing}"

    def test_output_schema(self, candidates):
        expected = {
            "rank", "target", "uniprot", "entrez_gene_symbol",
            "tier", "dual_somamer", "in_axon_pathway", "rationale",
        }
        assert expected.issubset(set(candidates.columns))

    def test_rank_is_sequential(self, candidates):
        assert list(candidates["rank"]) == list(range(1, len(candidates) + 1))

    def test_no_duplicate_targets(self, candidates):
        assert candidates["target"].is_unique

    def test_pathway_gene_extraction(self):
        if not ENRICHMENT_DIR.exists():
            pytest.skip("Step 2 output not found")
        genes = mod01.extract_pathway_genes(ENRICHMENT_DIR)
        assert len(genes) > 0
        # Known axon guidance genes should be present
        assert "PLXNA1" in genes or "SEMA3A" in genes or "NRP1" in genes


# --- Test sequence handling ---

class TestSequences:

    def test_parse_fasta_sequence(self):
        fasta = ">sp|P15692|VEGFA_HUMAN\nMNFLLSWVH\nWSLALLYLH\n"
        seq = mod02.parse_fasta_sequence(fasta)
        assert seq == "MNFLLSWVHWSLALLYLH"

    def test_trim_sequence(self):
        seq = "ABCDEFGHIJ"
        assert mod02.trim_sequence(seq, 1, 5) == "ABCDE"
        assert mod02.trim_sequence(seq, 3, 7) == "CDEFG"

    def test_two_chain_fasta_format(self):
        with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as f:
            path = Path(f.name)

        sflt1_seq = "MVSYW" * 20  # 100 aa mock
        partner_seq = "ACDEF" * 10  # 50 aa mock

        mod02.write_two_chain_fasta(
            sflt1_seq, partner_seq, "TestProtein", "Q12345", path
        )

        content = path.read_text()
        lines = content.strip().split("\n")

        # Should have exactly 2 header lines
        headers = [l for l in lines if l.startswith(">")]
        assert len(headers) == 2, f"Expected 2 headers, got {len(headers)}"

        # First chain should be sFLT1
        assert "sFLT1" in headers[0]
        assert "P17948" in headers[0]

        # Second chain should be partner
        assert "TestProtein" in headers[1]
        assert "Q12345" in headers[1]

        # Sequences should be present
        seq_lines = [l for l in lines if not l.startswith(">")]
        total_seq = "".join(seq_lines)
        assert len(total_seq) == 150  # 100 + 50

        path.unlink()


# --- Test LSF script generation ---


def _make_mock_jobs(tmpdir):
    """Create mock candidates + FASTA and generate LSF scripts."""
    cands = pd.DataFrame({
        "target": ["VEGFA", "NRP1"],
        "uniprot": ["P15692", "O14786"],
        "entrez_gene_symbol": ["VEGFA", "NRP1"],
        "tier": [0, 1],
        "dual_somamer": [False, True],
        "in_axon_pathway": [True, True],
        "rationale": ["positive control", "discovery"],
        "rank": [1, 2],
    })
    cands_path = tmpdir / "candidates.csv"
    cands.to_csv(cands_path, index=False)

    fasta_dir = tmpdir / "fasta"
    fasta_dir.mkdir()
    (fasta_dir / "sflt1_vs_VEGFA.fasta").write_text(">A\nMVSYW\n>B\nACDEF\n")
    (fasta_dir / "sflt1_vs_NRP1.fasta").write_text(">A\nMVSYW\n>B\nGHIJK\n")

    jobs_dir = tmpdir / "jobs"
    scripts = mod03.generate_lsf_scripts(cands_path, fasta_dir, jobs_dir)
    return scripts, jobs_dir


class TestLSFScripts:

    def test_standalone_lsf_has_required_directives(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            scripts, jobs_dir = _make_mock_jobs(tmpdir)

            assert len(scripts) == 2
            content = scripts[0].read_text()

            # Required LSF directives
            assert "#BSUB -J af2_VEGFA" in content
            assert "#BSUB -q gpu" in content
            assert "#BSUB -n 4" in content
            assert "ngpus_excl_p=1" in content
            assert "#BSUB -W 24:00" in content

    def test_standalone_lsf_uses_wrapper(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            scripts, jobs_dir = _make_mock_jobs(tmpdir)
            content = scripts[0].read_text()

            # Should use wrapper-based execution, not inline AF command
            assert "AF_JOB_COMMAND_B64" in content
            assert "AF_JOB_NAME" in content
            assert "AF_SENTINEL_DIR" in content
            assert "wrapper.sh" in content

    def test_manifest_json_generated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _make_mock_jobs(tmpdir)
            manifest_path = tmpdir / "jobs" / "manifest.json"
            assert manifest_path.exists()

            manifest = json.loads(manifest_path.read_text())
            assert "af2_VEGFA" in manifest
            assert "af2_NRP1" in manifest
            assert manifest["af2_VEGFA"]["target"] == "VEGFA"
            assert manifest["af2_VEGFA"]["uniprot"] == "P15692"
            assert "command_b64" in manifest["af2_VEGFA"]

    def test_manifest_command_b64_decodes(self):
        import base64
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _make_mock_jobs(tmpdir)
            manifest_path = tmpdir / "jobs" / "manifest.json"
            manifest = json.loads(manifest_path.read_text())

            cmd_b64 = manifest["af2_VEGFA"]["command_b64"]
            decoded = base64.b64decode(cmd_b64).decode("utf-8")
            assert "run_alphafold.py" in decoded
            assert "multimer" in decoded
            assert "sflt1_vs_VEGFA.fasta" in decoded

    def test_wrapper_script_generated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _make_mock_jobs(tmpdir)
            wrapper = tmpdir / "jobs" / "wrapper.sh"
            assert wrapper.exists()

            content = wrapper.read_text()
            assert "AF_JOB_COMMAND_B64" in content
            assert "AF_SENTINEL_DIR" in content
            assert "completed.log" in content
            # EXIT trap for sentinel
            assert "trap" in content

    def test_orchestrator_lsf_generated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _make_mock_jobs(tmpdir)
            orch = tmpdir / "jobs" / "orchestrator.lsf"
            assert orch.exists()

            content = orch.read_text()
            assert "#BSUB -J af2_orchestrator" in content
            assert "barrier_wait" in content
            assert "submit_batch" in content
            assert "check_upstream_failures" in content
            assert "manifest_job_tsv" in content
            assert "af2_VEGFA" in content
            assert "af2_NRP1" in content

    def test_submit_orchestrator_generated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            _make_mock_jobs(tmpdir)
            submit = tmpdir / "jobs" / "submit_orchestrator.sh"
            assert submit.exists()

            content = submit.read_text()
            assert "--dry-run" in content
            assert "bsub" in content

    def test_skips_candidate_without_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            cands = pd.DataFrame({
                "target": ["VEGFA", "MISSING"],
                "uniprot": ["P15692", "XXXXX"],
                "entrez_gene_symbol": ["VEGFA", "MISSING"],
                "tier": [0, 1],
                "dual_somamer": [False, True],
                "in_axon_pathway": [True, True],
                "rationale": ["ctrl", "disc"],
                "rank": [1, 2],
            })
            cands_path = tmpdir / "candidates.csv"
            cands.to_csv(cands_path, index=False)

            fasta_dir = tmpdir / "fasta"
            fasta_dir.mkdir()
            (fasta_dir / "sflt1_vs_VEGFA.fasta").write_text(">A\nMVSYW\n>B\nACDEF\n")
            # No FASTA for MISSING

            jobs_dir = tmpdir / "jobs"
            scripts = mod03.generate_lsf_scripts(cands_path, fasta_dir, jobs_dir)

            assert len(scripts) == 1
            assert scripts[0].name == "af2_VEGFA.lsf"

            manifest = json.loads((jobs_dir / "manifest.json").read_text())
            assert "af2_VEGFA" in manifest
            assert "af2_MISSING" not in manifest


# --- Test result parsing ---

class TestParseResults:

    def test_iptm_extraction_from_mock(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create mock ranking_debug.json
            result_dir = tmpdir / "VEGFA"
            result_dir.mkdir()

            mock_ranking = {
                "order": ["model_1_multimer_v3_pred_0", "model_2_multimer_v3_pred_0"],
                "iptm+ptm": {
                    "model_1_multimer_v3_pred_0": 0.85,
                    "model_2_multimer_v3_pred_0": 0.78,
                },
                "iptm": {
                    "model_1_multimer_v3_pred_0": 0.82,
                    "model_2_multimer_v3_pred_0": 0.75,
                },
            }
            with open(result_dir / "ranking_debug.json", "w") as f:
                json.dump(mock_ranking, f)

            scores = mod04.parse_ranking_debug(result_dir / "ranking_debug.json")
            assert len(scores) == 2
            assert scores["model_1_multimer_v3_pred_0"]["iptm"] == 0.85

    def test_interchain_pae_computation(self):
        # Mock PAE matrix: 10 residues chain A, 8 residues chain B
        n = 18
        pae = np.random.uniform(5, 25, (n, n))

        # Make intra-chain blocks low PAE
        pae[:10, :10] = 3.0
        pae[10:, 10:] = 4.0

        # Inter-chain blocks high PAE
        pae[:10, 10:] = 20.0
        pae[10:, :10] = 22.0

        mean_pae = mod04.compute_interchain_pae(pae, chain_a_len=10)

        # Should be approximately 21.0 (average of 20 and 22)
        assert 20.0 < mean_pae < 22.0

    def test_interchain_pae_edge_cases(self):
        pae = np.ones((10, 10))
        # chain_a_len == total length -> nan
        assert np.isnan(mod04.compute_interchain_pae(pae, chain_a_len=10))
        # chain_a_len == 0 -> nan
        assert np.isnan(mod04.compute_interchain_pae(pae, chain_a_len=0))

    def test_interaction_classification(self):
        # High confidence
        assert mod04.classify_interaction(0.85, 8.0) == "high_confidence"
        # Predicted
        assert mod04.classify_interaction(0.65, 12.0) == "predicted"
        # Low confidence (ipTM too low)
        assert mod04.classify_interaction(0.4, 20.0) == "low_confidence"
        # Low confidence (PAE too high)
        assert mod04.classify_interaction(0.7, 18.0) == "low_confidence"
        # No data
        assert mod04.classify_interaction(np.nan, 10.0) == "no_data"

    def test_pae_json_loading(self):
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
            pae_data = [{"predicted_aligned_error": [[1.0, 2.0], [3.0, 4.0]]}]
            json.dump(pae_data, f)
            path = Path(f.name)

        pae = mod04.load_pae_matrix(path)
        assert pae is not None
        assert pae.shape == (2, 2)
        assert pae[0, 0] == 1.0
        assert pae[1, 1] == 4.0
        path.unlink()
