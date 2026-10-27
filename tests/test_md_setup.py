"""Synthetic MD setup checks; never invoke VMD, NAMD, or a real scheduler."""

from __future__ import annotations

import contextlib
import io
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from stage4_md_stability import configure, prepare, slurm
from stage4_md_stability.common import TEMPLATES, restart_step, tcl_word, write_json


def atom(serial: int, residue: int, x: float, *, name: str = "CA", chain: str = "A", insertion: str = "") -> str:
    return (f"ATOM  {serial:5d}  {name:<3} ALA {chain}{residue:4d}{insertion:1s}   "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}           C\n")


def fixture(directory: Path, *, hmr: bool = False) -> tuple[Path, Path]:
    case_dir = directory / "candidate_001"
    case_dir.mkdir()
    for folder in ("conf", "out", "logs"):
        (case_dir / folder).mkdir()
    for name in ("system.psf", "system.pdb"):
        (case_dir / name).write_text("synthetic setup placeholder\n")
    write_json(case_dir / "case.json", {
        "schema_version": 1, "state": "prepared", "sequence_id": "candidate_001",
        "box_lengths": [100.0, 100.0, 100.0], "box_origin": [0.0, 0.0, 0.0],
        "hmr": hmr, "hydrogen_mass": 3.024 if hmr else None,
    })
    manifest = directory / "systems.json"
    write_json(manifest, {"schema_version": 1, "systems": [{"sequence_id": "candidate_001"}]})
    return manifest, case_dir


def restart(prefix: Path, step: int = 10000) -> None:
    Path(str(prefix) + ".coor").write_text("synthetic coordinates\n")
    Path(str(prefix) + ".vel").write_text("synthetic velocities\n")
    Path(str(prefix) + ".xsc").write_text(f"# NAMD extended system\n# columns\n{step} 100 0 0 0 100 0 0 0 100 0 0 0\n")


class QuietTest(unittest.TestCase):
    def invoke(self, main, arguments):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            result = main([str(value) for value in arguments])
        return result, out.getvalue(), err.getvalue()


class PreparationTests(QuietTest):
    def test_residue_map_preserves_source_identity_and_selects_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / "source.pdb"
            source.write_text(atom(1, 10, 0.0) + atom(2, 10, 1.0, insertion="A")
                              + atom(3, 50, 2.0, chain="B"))
            output = directory / "normalized.pdb"
            rows = prepare.normalize_pdb(source, output, "A", "candidate_001")
            self.assertEqual([row["prepared_residue_number"] for row in rows], [1, 2])
            self.assertEqual([row["input_insertion_code"] for row in rows], ["", "A"])
            self.assertEqual([line[22:26].strip() for line in output.read_text().splitlines()
                              if line.startswith("ATOM")], ["1", "2"])
            self.assertNotIn(" B  50", output.read_text())

    def test_multiple_models_and_disulfide_records_are_rejected(self):
        for content in ("MODEL        1\n" + atom(1, 1, 0.0) + "ENDMDL\nMODEL        2\n",
                        "SSBOND synthetic record\n" + atom(1, 1, 0.0)):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp) / "input.pdb"
                source.write_text(content)
                with self.assertRaisesRegex(ValueError, "multiple models|disulfide"):
                    prepare.normalize_pdb(source, Path(tmp) / "protein.pdb", "A", "candidate_001")

    def test_mmcif_and_nonfinite_coordinates_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / "model.cif"
            source.write_text("synthetic placeholder")
            with self.assertRaisesRegex(ValueError, "requires PDB"):
                prepare.normalize_pdb(source, directory / "protein.pdb", "A", "candidate_001")
            source = directory / "model.pdb"
            source.write_text(atom(1, 1, float("nan")))
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                prepare.normalize_pdb(source, directory / "protein.pdb", "A", "candidate_001")

    def test_chain_breaks_are_not_joined_into_a_single_segment(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / "input.pdb"
            source.write_text(atom(1, 1, 0.0) + "TER\n" + atom(2, 2, 2.0))
            with self.assertRaisesRegex(ValueError, "multiple fragments"):
                prepare.normalize_pdb(source, directory / "protein.pdb", "A", "candidate_001")

    def test_selecting_a_later_chain_does_not_inherit_an_earlier_ter(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / "input.pdb"
            source.write_text(atom(1, 1, 0.0, chain="A") + "TER\n" + atom(2, 4, 2.0, chain="B"))
            rows = prepare.normalize_pdb(source, directory / "protein.pdb", "B", "candidate_001")
            self.assertEqual([row["input_residue_number"] for row in rows], ["4"])

    def test_preparation_records_success_only_after_vmd_outputs(self):
        with tempfile.TemporaryDirectory(prefix="md setup with spaces ") as tmp:
            directory = Path(tmp)
            source, topology = directory / "input.pdb", directory / "topology.rtf"
            source.write_text(atom(1, 30, 0.0))
            topology.write_text("synthetic topology placeholder\n")
            output = directory / "results"

            def fake_vmd(command, *, cwd, stdout, stderr):
                self.assertEqual(command[1:], ["-dispdev", "text", "-e", "build.tcl"])
                for name in ("system.psf", "system.pdb"):
                    (cwd / name).write_text("synthetic placeholder\n")
                (cwd / "box.dat").write_text("100 100 100 0 0 0\n")
                stdout.write("STAGE4_PREPARED\n")
                return subprocess.CompletedProcess(command, 0)

            arguments = [source, "--source", "af3", "--sequence-id", "candidate_001",
                         "--topology", topology, "--outdir", output, "--hmr"]
            with patch.object(prepare.shutil, "which", return_value="mock-vmd"), \
                    patch.object(prepare.subprocess, "run", side_effect=fake_vmd):
                result, _, err = self.invoke(prepare.main, arguments)
                self.assertEqual(result, 0, err)
                result, _, err = self.invoke(prepare.main, arguments)
                self.assertEqual(result, 2)
                self.assertIn("new or empty", err)
            metadata = json.loads((output / "candidate_001/case.json").read_text())
            self.assertTrue(metadata["hmr"])
            self.assertEqual(metadata["structure_source"], "af3")
            self.assertEqual(metadata["state"], "prepared")
            self.assertTrue((output / "systems.json").exists())
            self.assertTrue((output / "candidate_001/out").is_dir())

    def test_vmd_exit_zero_without_success_marker_is_not_a_prepared_system(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source, topology = directory / "input.pdb", directory / "topology.rtf"
            source.write_text(atom(1, 1, 0.0))
            topology.write_text("synthetic placeholder\n")
            output = directory / "results"
            with patch.object(prepare.shutil, "which", return_value="mock-vmd"), \
                    patch.object(prepare.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)):
                result, _, err = self.invoke(prepare.main, [source, "--source", "design",
                    "--topology", topology, "--outdir", output])
            self.assertEqual(result, 2)
            self.assertIn("preparation failed", err)
            self.assertFalse((output / "systems.json").exists())
            self.assertFalse((output / "input/case.json").exists())

    def test_case_csv_resolves_paths_and_rejects_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "input.pdb").write_text(atom(1, 1, 0.0))
            cases = directory / "cases.csv"
            cases.write_text("sequence_id,input_structure,structure_source\ncandidate_001,input.pdb,af3\ncandidate_001,input.pdb,design\n")
            with patch.object(prepare.shutil, "which", return_value="mock-vmd"):
                result, _, err = self.invoke(prepare.main, ["--case-csv", cases,
                    "--topology", directory / "unused.rtf", "--outdir", directory / "output"])
            self.assertEqual(result, 2)
            self.assertIn("duplicate sequence_id", err)
            self.assertFalse((directory / "output").exists())


class ConfigurationTests(QuietTest):
    def test_initial_segment_uses_two_fs_and_creates_output_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            manifest, case = fixture(directory)
            shutil.rmtree(case / "out")
            result, _, err = self.invoke(configure.main, ["--manifest", manifest, "--segment", "equilibration",
                "--run-steps", "250000", "--parameter", directory / "parameter files/protein.prm"])
            self.assertEqual(result, 0, err)
            text = (case / "conf/equilibration.conf").read_text()
            self.assertIn("set timestep_fs 2.0", text)
            self.assertIn("Additional MD time: 0.5 ns", text)
            self.assertIn("interval: 50 ps", text)
            self.assertTrue((case / "out").is_dir())

    def test_four_fs_requires_hmr(self):
        for hmr, expected in ((False, 2), (True, 0)):
            with self.subTest(hmr=hmr), tempfile.TemporaryDirectory() as tmp:
                manifest, _ = fixture(Path(tmp), hmr=hmr)
                result, _, _ = self.invoke(configure.main, ["--manifest", manifest, "--segment", "initial",
                    "--run-steps", "10000", "--timestep-fs", "4", "--parameter", "/external/protein.prm"])
                self.assertEqual(result, expected)

    def test_restart_validates_all_files_and_uses_new_segment(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, case = fixture(Path(tmp))
            arguments = ["--manifest", manifest, "--segment", "production", "--restart-from", "equilibration",
                         "--run-steps", "10000", "--parameter", "/external/protein.prm"]
            result, _, err = self.invoke(configure.main, arguments)
            self.assertEqual(result, 2)
            self.assertIn("missing or empty", err)
            self.assertFalse((case / "conf/production.conf").exists())
            restart(case / "out/equilibration.restart", step=12340)
            result, _, err = self.invoke(configure.main, arguments)
            self.assertEqual(result, 0, err)
            text = (case / "conf/production.conf").read_text()
            self.assertIn("set minimize_steps 0", text)
            self.assertIn("Restart step at configuration time: 12340", text)
            self.assertIn('set restart_prefix "out/equilibration.restart"', text)

    def test_existing_configuration_or_outputs_are_not_overwritten(self):
        for existing in ("conf/initial.conf", "out/initial.dcd"):
            with self.subTest(existing=existing), tempfile.TemporaryDirectory() as tmp:
                manifest, case = fixture(Path(tmp))
                (case / existing).write_text("preserve me")
                result, _, _ = self.invoke(configure.main, ["--manifest", manifest, "--segment", "initial",
                    "--run-steps", "10000", "--parameter", "/external/protein.prm"])
                self.assertEqual(result, 2)
                self.assertEqual((case / existing).read_text(), "preserve me")

    def test_bad_cycle_counts_and_restart_minimization_are_rejected(self):
        for extra in (["--run-steps", "11"], ["--run-steps", "100", "--restart-from", "previous", "--minimize-steps", "200"]):
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as tmp:
                manifest, _ = fixture(Path(tmp))
                result, _, _ = self.invoke(configure.main, ["--manifest", manifest, "--segment", "initial",
                    "--parameter", "/external/protein.prm", *extra])
                self.assertEqual(result, 2)

    def test_malformed_restart_cell_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "previous.restart"
            restart(prefix)
            Path(str(prefix) + ".xsc").write_text("not valid XSC\n")
            with self.assertRaisesRegex(ValueError, "invalid NAMD restart"):
                restart_step(prefix)


@unittest.skipUnless(shutil.which("tclsh"), "Tcl helper checks require tclsh")
class TclTests(QuietTest):
    def tcl(self, code):
        # Tcl can print script errors while exiting zero; wrap the entire check.
        wrapped = "if {[catch {\n" + code + "\n} message]} {puts stderr $message; exit 1}\n"
        return subprocess.run(["tclsh"], input=wrapped, capture_output=True, text=True, check=False)

    def test_paths_roundtrip_without_tcl_substitution(self):
        value = 'directory with spaces/[not_a_command]/$not_a_variable/"quotes"/back\\slash'
        result = self.tcl(f"puts -nonewline {tcl_word(value)}")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, value)

    def test_asymmetric_coordinates_keep_full_padding(self):
        code = f"source {tcl_word(TEMPLATES / 'geometry.tcl')}\n"
        code += "lassign [md_water_box {0 0 0} {100 4 8} 25] center lengths lower upper\n"
        code += 'puts "$center|$lengths|$lower|$upper"'
        result = self.tcl(code)
        self.assertEqual(result.returncode, 0, result.stderr)
        center, lengths, lower, upper = [[float(v) for v in part.split()] for part in result.stdout.strip().split("|")]
        self.assertEqual(lengths, [150, 54, 58])
        for i, high in enumerate((100, 4, 8)):
            self.assertEqual((-center[i]) - lower[i], 25)
            self.assertEqual(upper[i] - (high - center[i]), 25)

    def test_hmr_checks_inner_bond_list_and_conserves_pair_mass(self):
        code = f"source {tcl_word(TEMPLATES / 'hmr.tcl')}\n"
        code += 'foreach neighbors {{} {7 8}} {if {![catch {md_single_partner [list $neighbors]}]} {error "invalid bond count accepted"}}\n'
        code += 'if {[md_single_partner {{7}}] != 7} {error "valid bond rejected"}\n'
        code += 'lassign [md_hmr_pair 1.008 12.011 3.024] h heavy\n'
        code += 'if {abs($h + $heavy - 13.019) > 1e-9} {error "mass not conserved"}\nputs $h'
        result = self.tcl(code)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(float(result.stdout.strip()), 3.024)

    def test_generated_namd_initial_and_restart_settings_execute_with_stubs(self):
        with tempfile.TemporaryDirectory(prefix="md paths [literal] ") as tmp:
            directory = Path(tmp)
            manifest, case = fixture(directory)
            parameter = directory / "parameters [literal] $value.prm"
            parameter.write_text("synthetic parameters\n")
            for segment, extra in (("initial", []), ("continuation", ["--restart-from", "initial"])):
                with self.subTest(segment=segment):
                    if extra:
                        restart(case / "out/initial.restart", step=12340)
                    result, _, err = self.invoke(configure.main, ["--manifest", manifest, "--segment", segment,
                        "--run-steps", "10000", "--parameter", parameter, *extra])
                    self.assertEqual(result, 0, err)
                    code = 'set calls {}\nproc unknown {args} {global calls; lappend calls $args}\n'
                    code += f"source {tcl_word(case / 'conf' / (segment + '.conf'))}\n"
                    code += 'foreach call $calls {puts [lindex $call 0]}\n'
                    result = self.tcl(code)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    commands = result.stdout.splitlines()
                    self.assertIn("run", commands)
                    if extra:
                        self.assertIn("binCoordinates", commands)
                        self.assertIn("binVelocities", commands)
                        self.assertIn("extendedSystem", commands)
                        self.assertNotIn("temperature", commands)
                        self.assertNotIn("minimize", commands)
                        self.assertNotIn("reinitvels", commands)
                    else:
                        self.assertIn("temperature", commands)
                        self.assertIn("minimize", commands)
                        self.assertIn("reinitvels", commands)

    def test_runtime_does_not_overwrite_existing_output_in_special_character_path(self):
        with tempfile.TemporaryDirectory(prefix="md [literal] ") as tmp:
            directory = Path(tmp)
            manifest, case = fixture(directory)
            parameter = directory / "parameters.prm"
            parameter.write_text("synthetic parameters\n")
            result, _, err = self.invoke(configure.main, ["--manifest", manifest, "--segment", "initial",
                "--run-steps", "10000", "--parameter", parameter])
            self.assertEqual(result, 0, err)
            (case / "out/initial.dcd").write_text("preserve me")
            result = self.tcl('proc unknown {args} {return}\n' +
                              f"source {tcl_word(case / 'conf/initial.conf')}")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("output prefix already has files", result.stderr)
            self.assertEqual((case / "out/initial.dcd").read_text(), "preserve me")


class SlurmTests(QuietTest):
    def test_default_only_writes_script_with_configurable_resources(self):
        with tempfile.TemporaryDirectory(prefix="md scheduler paths ") as tmp:
            directory = Path(tmp)
            _, case = fixture(directory)
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            with patch.object(slurm.subprocess, "run") as run:
                result, out, err = self.invoke(slurm.main, [config, "--namd", "/cluster path/namd3",
                    "--time", "24:00:00", "--cpus", "8", "--gpus", "2", "--account", "example",
                    "--partition", "gpu", "--gpu-type", "TEST", "--dependency", "123:456"])
                self.assertEqual(result, 0, err)
                run.assert_not_called()
            self.assertIn("Preview only", out)
            script = config.with_suffix(".slurm")
            text = script.read_text()
            self.assertIn("--nodes=1", text)
            self.assertIn("--ntasks-per-node=8", text)
            self.assertIn("--gres=gpu:TEST:2", text)
            self.assertIn("--dependency=afterok:123:456", text)
            self.assertIn("+p8", text)
            # Affinity and explicit devices are opt-in; see the two tests below.
            self.assertNotIn("+setcpuaffinity", text)
            self.assertNotIn("+devices", text)
            checked = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
            self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_cpu_affinity_and_devices_are_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, case = fixture(Path(tmp))
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            result, _, err = self.invoke(slurm.main, [config, "--time", "1:00:00", "--cpus", "8",
                "--gpus", "2", "--cpu-affinity", "--gpu-devices", "0,1"])
            self.assertEqual(result, 0, err)
            text = config.with_suffix(".slurm").read_text()
            self.assertIn("+setcpuaffinity", text)
            self.assertIn("+devices 0,1", text)

    def test_no_partition_is_emitted_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, case = fixture(Path(tmp))
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            result, _, err = self.invoke(slurm.main, [config, "--time", "1:00:00"])
            self.assertEqual(result, 0, err)
            text = config.with_suffix(".slurm").read_text()
            self.assertNotIn("--partition", text)
            self.assertIn("--qos=inferno", text)

    def test_submit_requires_explicit_flag_and_passes_arguments_without_shell(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, case = fixture(Path(tmp))
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            with patch.object(slurm.shutil, "which", return_value="mock-sbatch"), \
                    patch.object(slurm.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "12345\n", "")) as run:
                result, out, err = self.invoke(slurm.main, [config, "--namd", "namd3", "--time", "01:00:00", "--submit"])
                self.assertEqual(result, 0, err)
                self.assertIn("12345", out)
                self.assertEqual(run.call_args.args[0], ["sbatch", "--parsable", str(config.resolve().with_suffix(".slurm"))])
                self.assertNotIn("shell", run.call_args.kwargs)

    def test_site_option_injection_and_existing_scripts_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, case = fixture(Path(tmp))
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            arguments = [config, "--namd", "namd3", "--time", "01:00:00"]
            result, _, _ = self.invoke(slurm.main, [*arguments, "--partition", "gpu\n#SBATCH --bad"])
            self.assertEqual(result, 2)
            script = config.with_suffix(".slurm")
            script.write_text("preserve me")
            result, _, _ = self.invoke(slurm.main, arguments)
            self.assertEqual(result, 2)
            self.assertEqual(script.read_text(), "preserve me")

    def test_preview_can_be_submitted_explicitly_with_identical_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, case = fixture(Path(tmp))
            config = case / "conf/initial.conf"
            config.write_text("synthetic configuration\n")
            arguments = [config, "--namd", "namd3", "--time", "01:00:00"]
            with patch.object(slurm.shutil, "which", return_value="mock-sbatch"), \
                    patch.object(slurm.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "12345\n", "")) as run:
                result, _, err = self.invoke(slurm.main, arguments)
                self.assertEqual(result, 0, err)
                run.assert_not_called()
                result, out, err = self.invoke(slurm.main, [*arguments, "--submit"])
                self.assertEqual(result, 0, err)
                self.assertIn("Reusing", out)
                run.assert_called_once()


if __name__ == "__main__":
    unittest.main()


class NamdTemplateTests(unittest.TestCase):
    """Pin the template to the proven base.conf from the previous version of this workflow.

    That configuration has run to completion on the target cluster; configurations that
    drifted from it crashed. A change here should be a deliberate decision, not a cleanup.
    """

    # Keyword -> value as written in base.conf (Tcl variables where the template takes
    # the value from configure.py; their defaults match base.conf).
    PROVEN = {
        "wrapAll": "on", "wrapNearest": "on", "exclude": "scaled1-4", "1-4scaling": "1.0",
        "cutoff": "12.0", "switching": "on", "switchdist": "10.0", "pairlistdist": "14.0",
        "vdwForceSwitching": "on", "rigidBonds": "all", "nonbondedFreq": "1",
        "fullElectFrequency": "2", "stepspercycle": "10", "PME": "yes", "PMEGridSpacing": "1.0",
        "langevin": "on", "langevinDamping": "1.0", "useGroupPressure": "no",
        "useFlexibleCell": "no", "useConstantArea": "no", "langevinPiston": "on",
        "langevinPistonPeriod": "100.0", "langevinPistonDecay": "50.0",
    }

    def setUp(self):
        self.settings = {}
        for line in (TEMPLATES / "namd.tcl").read_text().splitlines():
            fields = line.split()
            if len(fields) == 2 and not line.startswith((" ", "#")):
                self.settings[fields[0]] = fields[1]

    def test_settings_match_the_proven_configuration(self):
        for keyword, value in self.PROVEN.items():
            self.assertEqual(self.settings.get(keyword), value, keyword)

    def test_no_settings_beyond_the_proven_configuration(self):
        # base.conf leaves these to NAMD's defaults; the template must not add them.
        for keyword in ("useSettle", "binaryoutput", "binaryrestart", "oneFourScaling", "margin"):
            self.assertNotIn(keyword, self.settings, keyword)

    def test_gpu_resident_is_opt_in(self):
        self.assertIn("if {$gpu_resident}", (TEMPLATES / "namd.tcl").read_text())

    def test_configure_defaults_match_the_proven_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, case = fixture(Path(tmp), hmr=True)
            self.assertEqual(configure.main(["--manifest", str(manifest), "--segment", "seg",
                                             "--run-steps", "1000"]), 0)
            text = (case / "conf/seg.conf").read_text()
        for line in ("set timestep_fs 4.0", "set minimize_steps 200", "set output_every 25000",
                     "set energy_every 5000", "set temperature_k 300.0",
                     "set pressure_bar 1.01325", "set gpu_resident 0"):
            self.assertIn(line, text)


class GpuResidentConfigureTests(unittest.TestCase):
    def test_flag_controls_the_emitted_setting(self):
        for flag, expected in ((None, "set gpu_resident 0"), ("--gpu-resident", "set gpu_resident 1")):
            with tempfile.TemporaryDirectory() as tmp:
                manifest, case = fixture(Path(tmp))
                argv = ["--manifest", str(manifest), "--segment", "seg", "--run-steps", "1000"]
                if flag:
                    argv.append(flag)
                self.assertEqual(configure.main(argv), 0)
                text = (case / "conf/seg.conf").read_text()
                self.assertIn(expected, text)
