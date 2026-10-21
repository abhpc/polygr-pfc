#!/usr/bin/env python3
"""Frozen CPU/DCU comparisons and small image checks; run only in a Debug job."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

import numpy as np
from PIL import Image

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "source"))
import plotter

BASELINE_COMMIT = "c86efb51cf03086fb98897ab8abfe2704a62096e"
DENSITY_ATOL, DENSITY_RTOL = 1e-6, 1e-5
SCIENCE_ATOL, SCIENCE_RTOL = 1e-8, 1e-6


def require(condition, message):
	if not condition:
		raise ValueError(message)


def sha256(path):
	digest = hashlib.sha256()
	with Path(path).open("rb") as stream:
		for chunk in iter(lambda: stream.read(1024 * 1024), b""):
			digest.update(chunk)
	return digest.hexdigest()


def error_metrics(reference, actual, atol, rtol):
	require(reference.shape == actual.shape, "comparison shape mismatch")
	require(np.isfinite(reference).all() and np.isfinite(actual).all(), "nonfinite comparison")
	difference = np.abs(actual - reference)
	ratio = difference / (atol + rtol * np.abs(reference))
	index = int(np.argmax(ratio))
	return {
		"max_absolute_error": float(difference.max()),
		"max_tolerance_ratio": float(ratio.max()),
		"over_limit": int(np.count_nonzero(ratio > 1)),
		"worst_flat_index": index,
		"reference_at_worst": float(reference.flat[index]),
		"actual_at_worst": float(actual.flat[index]),
	}


def read_science(path):
	data = np.loadtxt(path, dtype=np.float64, comments=None, ndmin=2)
	require(data.shape[1] == 6 and data.shape[0] > 0, f"{path}: expected six scientific columns")
	require(np.isfinite(data).all(), f"{path}: nonfinite scientific data")
	require(np.equal(data[:, :2], np.floor(data[:, :2])).all(), f"{path}: step/time must be integers")
	require((data[:, :2] >= 0).all(), f"{path}: negative step/time")
	return data


def input_configuration(path):
	tokens = []
	for line in Path(path).read_text().splitlines():
		tokens.extend(line.split("#", 1)[0].split())
	position, blocks, print_interval, write_interval = 0, [], None, None
	while position < len(tokens):
		label = tokens[position]
		position += 1
		if label in ("S", "O", "A", "M"):
			count = {"S": 1, "O": 2, "A": 2, "M": 5}[label]
			if label == "O":
				print_interval, write_interval = map(int, tokens[position:position + count])
			position += count
		elif label == "I":
			kind = int(tokens[position])
			position += 1 + {0: 2, 1: 7, 2: 3}[kind]
		elif label == "R":
			values = tokens[position:position + 5]
			require(len(values) == 5 and print_interval and write_interval, f"{path}: incomplete R/O input")
			blocks.append({"steps": int(values[0]), "optimize": int(values[4]),
				"print": print_interval, "write": write_interval})
			position += 5
		else:
			raise ValueError(f"{path}: unknown input label {label}")
	return blocks


def records(path, prefix):
	result = []
	for line in Path(path).read_text(errors="replace").splitlines():
		if line.startswith(prefix + " "):
			fields = dict(token.split("=", 1) for token in shlex.split(line)[1:])
			result.append(fields)
	return result


def check_timing(path, phase, blocks):
	entries = records(path, "PFC_TIMING")
	require(entries, f"{path}: missing PFC_TIMING")
	for entry in entries:
		require(entry.get("phase") == phase, f"{path}: timing phase mismatch")
		seconds = float(entry["seconds"])
		require(math.isfinite(seconds) and seconds >= 0, f"{path}: invalid timing seconds")
		require(entry["kind"] in {"init", "step", "optimize", "total"}, f"{path}: unknown timing kind")
	for kind in ("init", "total"):
		require(sum(entry["kind"] == kind for entry in entries) == 1, f"{path}: expected one {kind} timing")
	require(sum(entry["kind"] == "step" for entry in entries) == sum(block["steps"] for block in blocks),
		f"{path}: extra/missing advancement timings")
	optimization_count = 0
	for block_index, block in enumerate(blocks, 1):
		steps = [entry for entry in entries if entry["kind"] == "step" and int(entry["block"]) == block_index]
		require([int(entry["step"]) for entry in steps] == list(range(1, block["steps"] + 1)),
			f"{path}: missing, duplicate or unordered step timings in block {block_index}")
		cumulative = [float(entry["cumulative_seconds"]) for entry in steps]
		require(all(math.isfinite(value) and value >= 0 for value in cumulative), f"{path}: invalid cumulative timing")
		require(all(a <= b for a, b in zip(cumulative, cumulative[1:])), f"{path}: decreasing cumulative timing")
		require(all(float(entry["seconds"]) <= value + 1e-8 for entry, value in zip(steps, cumulative)),
			f"{path}: advancement duration exceeds cumulative time")
		interval = block["optimize"]
		expected = list(range(interval, block["steps"] + 1, interval)) if interval > 0 else []
		optimized = [entry for entry in entries if entry["kind"] == "optimize" and int(entry["block"]) == block_index]
		require([int(entry["step"]) for entry in optimized] == expected,
			f"{path}: optimization timing schedule differs from original input")
		optimization_count += len(expected)
	require(sum(entry["kind"] == "optimize" for entry in entries) == optimization_count,
		f"{path}: extra optimization timings")
	return {"advancement_records": sum(block["steps"] for block in blocks), "optimization_records": optimization_count,
		"init_seconds": float(next(entry["seconds"] for entry in entries if entry["kind"] == "init")),
		"total_seconds": float(next(entry["seconds"] for entry in entries if entry["kind"] == "total"))}


def slab_table(path):
	rows = []
	for line in Path(path).read_text().splitlines():
		if not line.strip() or line.lstrip().startswith("#"):
			continue
		values = line.split()
		if values[0] == "PFC_SLAB":
			entry = dict(value.split("=", 1) for value in values[1:])
			rows.append(tuple(int(entry[key]) for key in ("rank", "lH", "lh0", "lW", "lw0")))
		else:
			require(len(values) == 5, f"{path}: expected rank lH lh0 lW lw0")
			rows.append(tuple(map(int, values)))
	return sorted(rows)


def check_coverage(rows, width, height):
	require([row[0] for row in rows] == list(range(32)), "slab evidence must contain each of 32 ranks once")
	for length_column, start_column, total in ((1, 2, height), (3, 4, width // 2 + 1)):
		spans = sorted((row[start_column], row[length_column]) for row in rows if row[length_column] > 0)
		end = 0
		for start, length in spans:
			require(start == end and length >= 0, "slab overlap/gap")
			end += length
		require(end == total, "slabs do not cover global field")


def check_devices(path, cpu_slabs, width, height, multinode):
	devices, slabs = records(path, "PFC_DEVICE"), records(path, "PFC_SLAB")
	execution = records(path, "PFC_EXEC")
	transfers = records(path, "PFC_TRANSFER")
	require(sorted(int(entry["rank"]) for entry in devices) == list(range(32)), f"{path}: missing/duplicate device ranks")
	actual = sorted(tuple(int(entry[key]) for key in ("rank", "lH", "lh0", "lW", "lw0")) for entry in slabs)
	check_coverage(actual, width, height)
	require(actual == cpu_slabs, f"{path}: CPU/DCU scientific slab mismatch")
	nodes = {}
	for entry in devices:
		host = entry.get("host")
		require(host and int(entry["visible"]) == 4, f"{path}: each node must expose four allocated devices")
		require(entry["role"] in ("leader", "member"), f"{path}: invalid execution role")
		nodes.setdefault(host, []).append(entry)
	require(len(nodes) == (2 if multinode else 1), f"{path}: wrong node count")
	expected_per_node = 16 if multinode else 32
	expected_leaders = 8 if multinode else 4
	groups = []
	slabs_by_rank = {row[0]: row for row in actual}
	for host, entries in nodes.items():
		require(sorted(int(entry["local_rank"]) for entry in entries) == list(range(expected_per_node)),
			f"{path}: local rank coverage mismatch on {host}")
		cpus = [int(entry["cpu"]) for entry in entries]
		require(all(cpu >= 0 for cpu in cpus) and len(set(cpus)) == expected_per_node,
			f"{path}: scientific MPI ranks must occupy distinct CPUs on {host}")
		cards = {}
		for entry in entries:
			device, local_rank = int(entry["device"]), int(entry["local_rank"])
			require(device == local_rank % 4, f"{path}: incorrect local rank/device mapping")
			cards.setdefault(entry["pci"], []).append(entry)
		require(len(cards) == 4 and all(len(values) == expected_per_node // 4 for values in cards.values()),
			f"{path}: four distinct physical DCUs must actually participate on {host}")
		require(all(len({int(entry["device"]) for entry in values}) == 1 for values in cards.values()),
			f"{path}: device identity collision")
		require(sum(entry["role"] == "leader" for entry in entries) == 4, f"{path}: each node needs four execution leaders")
		for pci, members in cards.items():
			leaders = [entry for entry in members if entry["role"] == "leader"]
			require(len(leaders) == 1, f"{path}: physical card needs exactly one execution leader")
			leader_rank = int(leaders[0]["rank"])
			require(all(int(entry["leader_rank"]) == leader_rank for entry in members),
				f"{path}: member points to a different host/card execution leader")
			ranks = sorted(int(entry["rank"]) for entry in members)
			groups.append({"host": host, "pci": pci, "leader_rank": leader_rank, "member_ranks": ranks,
				"lH": sum(slabs_by_rank[rank][1] for rank in ranks),
				"lW": sum(slabs_by_rank[rank][3] for rank in ranks)})
	require(len(execution) == expected_leaders and len({int(entry["leader_rank"]) for entry in execution}) == expected_leaders,
		f"{path}: missing/duplicate PFC_EXEC records")
	execution_by_rank = {int(entry["leader_rank"]): entry for entry in execution}
	require(set(execution_by_rank) == {group["leader_rank"] for group in groups}, f"{path}: execution leader identities differ")
	require(len(transfers) == expected_leaders and len({int(entry["leader_rank"]) for entry in transfers}) == expected_leaders,
		f"{path}: missing/duplicate PFC_TRANSFER records")
	transfer_by_rank = {int(entry["leader_rank"]): entry["mode"] for entry in transfers}
	require(set(transfer_by_rank) == set(execution_by_rank), f"{path}: transfer leader identities differ")
	modes = set(transfer_by_rank.values())
	expected_mode = "host-mpi" if multinode else "host-shared"
	require(modes == {expected_mode}, f"{path}: expected {expected_mode} transfers")
	for group in groups:
		entry = execution_by_rank[group["leader_rank"]]
		group["transfer_mode"] = transfer_by_rank[group["leader_rank"]]
		require(int(entry["leaders"]) == expected_leaders and int(entry["members"]) == expected_per_node // 4,
			f"{path}: incorrect execution communicator/group size")
		require(int(entry["lH"]) == group["lH"] and int(entry["lW"]) == group["lW"],
			f"{path}: execution aggregate differs from original member slabs")
	require(sorted(rank for group in groups for rank in group["member_ranks"]) == list(range(32)),
		f"{path}: original scientific ranks must each belong to exactly one execution group")
	return {"hosts": {host: sorted({entry["pci"] for entry in entries}) for host, entries in nodes.items()},
		"execution_groups": groups,
		"cpu_affinity": {host: {str(entry["rank"]): int(entry["cpu"]) for entry in entries} for host, entries in nodes.items()}}


def check_metadata(cpu, dcu, width, height, multinode):
	metadata = {}
	for mode, directory in (("cpu", cpu), ("dcu", dcu)):
		path = directory.parent / "metadata.json"
		info = json.loads(path.read_text())
		for key in ("baseline_commit", "baseline_source_sha256", "cpu_binary_sha256", "dcu_binary_sha256",
			"host_flags", "fftw_version", "fftw_config", "mpi_size", "nodes", "partition", "job_id",
			"exit_status", "build_evidence", "fixture_hashes"):
			require(key in info, f"{path}: missing evidence key {key}")
		require(info["baseline_commit"] == BASELINE_COMMIT, f"{path}: original CPU commit mismatch")
		require(info["partition"] == "Debug" and int(info["mpi_size"]) == 32 and int(info["exit_status"]) == 0,
			f"{path}: invalid job resources/status")
		require(int(info["nodes"]) == (2 if mode == "dcu" and multinode else 1), f"{path}: wrong metadata node count")
		require(info["job_id"], f"{path}: missing Slurm job identity")
		require(info["fftw_version"] == "3.3.10", f"{path}: unverified FFTW version")
		for flag in ("-O3", "-march=native", "-ffp-contract=off", "-Wall"):
			require(flag in info["host_flags"], f"{path}: missing shared host flag {flag}")
		for flag in ("--enable-mpi", "--enable-sse2", "--enable-avx", "--enable-avx2"):
			require(flag in info["fftw_config"], f"{path}: missing FFTW optimization {flag}")
		for key in ("baseline_source_sha256", "cpu_binary_sha256", "dcu_binary_sha256"):
			require(re.fullmatch(r"[0-9a-f]{64}", info[key]) is not None, f"{path}: invalid {key}")
		build_path = Path(info["build_evidence"])
		require(build_path.is_file(), f"{path}: missing build evidence")
		build = json.loads(build_path.read_text())
		for key in ("baseline_commit", "baseline_source_sha256", "cpu_binary_sha256", "host_flags", "fftw_version", "fftw_config"):
			require(build[key] == info[key], f"{path}: recorded CPU reference differs from build evidence: {key}")
		if mode == "dcu":
			require(build["dcu_binary_sha256"] == info["dcu_binary_sha256"], f"{path}: DCU binary differs from build evidence")
		fixture = REPOSITORY / "tests" / "fixtures" / ("v3" if width == height == 4240 else "v1")
		for phase in ("step1", "step2"):
			name = phase + ".in"
			digest = sha256(directory / name)
			require(digest == sha256(fixture / name) == info["fixture_hashes"][name], f"{path}: frozen input hash mismatch {name}")
		metadata[mode] = info
	for key in ("baseline_source_sha256", "cpu_binary_sha256", "host_flags", "fftw_config"):
		require(metadata["cpu"][key] == metadata["dcu"][key], f"CPU/DCU build evidence differs: {key}")
	return metadata


def validate(args, evidence):
	cpu, dcu = args.cpu.resolve(), args.dcu.resolve()
	require((args.width, args.height) in ((256, 256), (4240, 4240)), "comparison requires a frozen V1/V2 or 100 nm V3 grid")
	if args.width == args.height == 4240:
		spacing_angstrom = .7 * 2.46 / 7.3
		evidence["physical_box"] = {"target_nm": [100, 100],
			"converter_scales": [.7, .7, 7.3, 2.46], "nominal_spacing_angstrom": spacing_angstrom,
			"nominal_size_nm": [args.width * spacing_angstrom / 10, args.height * spacing_angstrom / 10]}
	if args.speedup_min is not None:
		require(not args.multinode and not args.reuse_cpu and cpu.parent == dcu.parent,
			"performance pair must share one single-node Slurm job directory")
	else:
		require(args.multinode or args.reuse_cpu or cpu.parent == dcu.parent,
			"cross-job V1 correctness comparison requires explicit --reuse-cpu")
	evidence["metadata"] = check_metadata(cpu, dcu, args.width, args.height, args.multinode)
	evidence["reused_cpu"] = bool(args.reuse_cpu or args.multinode)
	if args.speedup_min is not None:
		require(evidence["metadata"]["cpu"]["job_id"] == evidence["metadata"]["dcu"]["job_id"],
			"performance CPU/DCU metadata must identify the same Slurm job")
		bootstrap = evidence["metadata"]["dcu"].get("mpi_bootstrap")
		require(bootstrap in ("fork", "slurm") and evidence["metadata"]["cpu"].get("mpi_bootstrap") == bootstrap,
			"performance CPU/DCU pair must use the same recorded MPI bootstrap")
	cpu_slabs = slab_table(cpu / "slabs.tsv")
	check_coverage(cpu_slabs, args.width, args.height)
	evidence["slabs"] = cpu_slabs
	evidence["phases"] = {}
	all_numeric_passed = True
	initialized_fields_equal = True
	expected_names = set()
	for phase in ("step1", "step2"):
		blocks = input_configuration(dcu / (phase + ".in"))
		require(len(blocks) == 1 and blocks[0]["steps"] == args.steps, "frozen case workload mismatch")
		expected_names.update(f"{phase}-{step}.dat" for block in blocks
			for step in range(0, block["steps"] + 1, block["write"]))
		result = {"devices": check_devices(dcu / (phase + ".stdout"), cpu_slabs,
			args.width, args.height, args.multinode), "timing": check_timing(dcu / (phase + ".stdout"), phase, blocks)}
		evidence["phases"][phase] = result
		for step, label in ((0, "initial"), (args.steps, "final")):
			reference_path, actual_path = cpu / f"{phase}-t-{step}.dat", dcu / f"{phase}-{step}.dat"
			reference = plotter.read_density(reference_path, args.width, args.height)
			actual = plotter.read_density(actual_path, args.width, args.height)
			metrics = error_metrics(reference, actual, DENSITY_ATOL, DENSITY_RTOL)
			metrics.update({"cpu_path": str(reference_path), "dcu_path": str(actual_path),
				"cpu_sha256": sha256(reference_path), "dcu_sha256": sha256(actual_path)})
			if phase == "step1" and step == 0:
				initialized_fields_equal = bool(np.array_equal(reference, actual))
				metrics["exactly_equal"] = initialized_fields_equal
				metrics["exact_bytes_equal"] = metrics["cpu_sha256"] == metrics["dcu_sha256"]
				initialized_fields_equal &= metrics["exact_bytes_equal"]
			result[label] = metrics
			all_numeric_passed &= metrics["over_limit"] == 0
		reference, actual = read_science(cpu / (phase + ".out")), read_science(dcu / (phase + ".out"))
		expected_steps = [step for block in blocks for step in range(0, block["steps"] + 1, block["print"])]
		require(reference[:, 0].tolist() == actual[:, 0].tolist() == expected_steps,
			f"{phase}: original scientific diagnostic schedule changed")
		result["science"] = {}
		for column, name in ((2, "dx"), (3, "dy"), (4, "free_energy"), (5, "mean_density")):
			metrics = error_metrics(reference[:, column], actual[:, column], SCIENCE_ATOL, SCIENCE_RTOL)
			result["science"][name] = metrics
			all_numeric_passed &= metrics["over_limit"] == 0
		result["cpu_seconds"] = float((cpu / (phase + ".time")).read_text().strip())
		result["dcu_seconds"] = float((dcu / (phase + ".time")).read_text().strip())
		require(all(math.isfinite(result[key]) and result[key] > 0 for key in ("cpu_seconds", "dcu_seconds")),
			f"{phase}: invalid complete external wall time")
		evidence["phases"][phase] = result
	require({path.name for path in dcu.glob("*.dat")} == expected_names, "unexpected or missing DCU density output names")
	log = (dcu / "step2.stdout").read_text()
	require("PFC_RESTART_PARSE mode=parallel" in log, "canonical restart did not use distributed parsing")
	require("step1-t-" + str(args.steps) + ".dat" in log and "step1-" + str(args.steps) + ".dat" in log,
		"missing documented old-name restart fallback record")
	counts = {}
	for mode, directory in (("cpu", cpu), ("dcu", dcu)):
		count = int((directory / "atoms.count").read_text().strip())
		require(count > 0, f"{directory}: final conversion yielded no atoms")
		suffix = f"step2-t-{args.steps}.xy" if mode == "cpu" else f"step2-{args.steps}.xy"
		coordinates = np.loadtxt(directory / suffix, ndmin=2)
		require(np.isfinite(coordinates).all() and coordinates.shape[0] == count,
			f"{directory}: converted coordinate count/nonfinite mismatch")
		lammps = (directory / "graphene.data").read_text()
		match = re.search(r"^\s*(\d+)\s+atoms\s*$", lammps, flags=re.MULTILINE)
		require(match and int(match.group(1)) == count, f"{directory}: LAMMPS atom count mismatch")
		counts[mode] = count
	evidence["atom_counts"] = counts
	require(counts["cpu"] == counts["dcu"], "CPU/DCU final atom counts differ")
	evidence["numerical_passed"] = bool(all_numeric_passed and initialized_fields_equal)
	require(initialized_fields_equal, "phase-one CPU/DCU original initialized fields differ")
	require(all_numeric_passed, "CPU/DCU field or diagnostic comparison exceeded frozen tolerance")
	cpu_seconds = sum(phase["cpu_seconds"] for phase in evidence["phases"].values())
	dcu_seconds = sum(phase["dcu_seconds"] for phase in evidence["phases"].values())
	evidence["performance"] = {"cpu_complete_seconds": cpu_seconds, "dcu_complete_seconds": dcu_seconds,
		"speedup": cpu_seconds / dcu_seconds, "required_speedup": args.speedup_min}
	if args.speedup_min is not None:
		require(cpu_seconds / dcu_seconds >= args.speedup_min, "complete two-stage speedup is below frozen requirement")


def images(args, evidence):
	root = args.root.resolve()
	root.mkdir(parents=True, exist_ok=False)
	tool = REPOSITORY / "source" / "plotter.py"
	checks = []

	def invoke(arguments, success=True, contains=None):
		result = subprocess.run([sys.executable, str(tool), *map(str, arguments)], text=True,
			stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		checks.append({"arguments": list(map(str, arguments)), "exit_status": result.returncode,
			"stdout": result.stdout, "stderr": result.stderr})
		require((result.returncode == 0) == success, f"plotter unexpected exit: {arguments}: {result.stderr}")
		if contains:
			require(contains in result.stdout + result.stderr, f"plotter missing expected message {contains}")
		return result

	def check_image(source, output, width, height):
		field = plotter.read_density(source, width, height)
		low, high = float(field.min()), float(field.max())
		expected = np.zeros((height, width), dtype=np.uint8) if high == low else \
			np.flipud(np.floor(255 * (field.reshape(height, width) - low) / (high - low) + .5)).astype(np.uint8)
		with Image.open(output) as image:
			require(image.mode == "L" and image.size == (width, height), f"{output}: invalid grayscale mode/dimensions")
			require(np.max(np.abs(np.asarray(image).astype(np.int16) - expected.astype(np.int16))) <= 1,
				f"{output}: grayscale/orientation differs from original convention")

	invoke([args.data / "step1-#.dat", root / "step1-#.png", args.width, args.height, 0, args.steps, args.steps],
		contains="PFC_PLOT generated=2 missing=0")
	for step in (0, args.steps):
		check_image(args.data / f"step1-{step}.dat", root / f"step1-{step}.png", args.width, args.height)
	invoke([args.data / f"step2-{args.steps}.dat", root / f"step2-{args.steps}.png", args.width, args.height])
	check_image(args.data / f"step2-{args.steps}.dat", root / f"step2-{args.steps}.png", args.width, args.height)
	fixture = root / "tiny"
	fixture.mkdir()
	(fixture / "a-0.dat").write_text("0\n1\n2\n3\n")
	(fixture / "constant-0.dat").write_text("7\n7\n7\n7\n")
	unsafe = "a:\u4e2d-0.png"
	safe = "a____-0.png"
	invoke([fixture / "a-t-0.dat", fixture / unsafe, 2, 2], contains="PFC_INPUT")
	check_image(fixture / "a-0.dat", fixture / safe, 2, 2)
	before = sha256(fixture / safe)
	invoke([fixture / "a-0.dat", fixture / unsafe, 2, 2], success=False)
	require(sha256(fixture / safe) == before, "existing PNG was overwritten")
	(fixture / "a-t-0.dat").write_text("3\n2\n1\n0\n")
	invoke([fixture / "a-t-0.dat", fixture / "exact-0.png", 2, 2])
	check_image(fixture / "a-t-0.dat", fixture / "exact-0.png", 2, 2)
	invoke([fixture / "constant-0.dat", fixture / "constant-0.png", 2, 2])
	with Image.open(fixture / "constant-0.png") as image:
		require(not np.asarray(image).any(), "constant field must be black")
	invoke([fixture / "a-#.dat", fixture / "missing-#.png", 2, 2, 0, 1, 1],
		contains="PFC_PLOT generated=1 missing=1")
	invoke([fixture / "absent-#.dat", fixture / "absent-#.png", 2, 2, 0, 1, 1], success=False)
	invoke([fixture / "absent-0.dat", fixture / "absent-0.png", 2, 2], success=False)
	for start, increment, end in ((0, 0, 1), (2, 1, 1), (-1, 1, 1)):
		invoke([fixture / "a-#.dat", fixture / "range-#.png", 2, 2, start, increment, end], success=False)
	for index, invalid in enumerate(("1\n2\n3\n", "1 2\n3 4\n", "1\n2\nnan\n4\n",
		"1\n2\ninf\n4\n", "2 2\n1\n2\n3\n4\n", "")):
		bad = fixture / f"bad{index}-0.dat"
		bad.write_text(invalid)
		invoke([bad, fixture / f"bad{index}-0.png", 2, 2], success=False)
		require(not (fixture / f"bad{index}-0.png").exists(), "corrupt field produced PNG")
	invoke([fixture / "a-0.dat", fixture / "wrong-1.png", 2, 2], success=False)
	invoke([fixture / "a-#.dat", fixture / "bad##-#.png", 2, 2, 0, 1, 1], success=False)
	for step in (2, 10):
		(fixture / f"order-{step}.dat").write_text("0\n1\n2\n3\n")
	result = invoke([fixture / "order-#.dat", fixture / "order-#.png", 2, 2, 2, 8, 10])
	ordered = [line for line in result.stdout.splitlines() if line.startswith("PFC_PLOT_FRAME ")]
	require(len(ordered) == 2 and "order-2.png" in ordered[0] and "order-10.png" in ordered[1],
		"batch images were not processed in integer step order")
	require(plotter.safe_prefix("x:\u4e2d") == "x____", "C/Python bytewise safety rule differs")
	for prefix in ("a:", "a?", "a\udcff"):
		require(plotter.safe_prefix(prefix) == "a_", "different invalid bytes must map to the same safe prefix")
	special_directory = fixture / "dir:\u4e2d"
	special_directory.mkdir()
	(special_directory / "field____-0.dat").write_text("0\n1\n2\n3\n")
	invoke([special_directory / "field:\u4e2d-t-0.dat", special_directory / "field:\u4e2d-0.png", 2, 2],
		contains="PFC_INPUT")
	check_image(special_directory / "field____-0.dat", special_directory / "field____-0.png", 2, 2)
	for invalid in ("",):
		try:
			plotter.safe_prefix(invalid)
		except ValueError:
			pass
		else:
			raise ValueError("empty safe prefix accepted")
	evidence["real_frames"] = 3
	evidence["checks"] = checks


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	subparsers = parser.add_subparsers(dest="command", required=True)
	comparison = subparsers.add_parser("validate")
	comparison.add_argument("--cpu", type=Path, required=True)
	comparison.add_argument("--dcu", type=Path, required=True)
	comparison.add_argument("--width", type=int, required=True)
	comparison.add_argument("--height", type=int, required=True)
	comparison.add_argument("--steps", type=int, required=True)
	comparison.add_argument("--multinode", action="store_true")
	comparison.add_argument("--reuse-cpu", action="store_true",
		help="reuse a traceable V1 CPU result after a targeted DCU fix; never for performance acceptance")
	comparison.add_argument("--speedup-min", type=float)
	comparison.add_argument("--evidence", type=Path, required=True)
	image_parser = subparsers.add_parser("images")
	image_parser.add_argument("--root", type=Path, required=True)
	image_parser.add_argument("--data", type=Path, required=True)
	image_parser.add_argument("--width", type=int, required=True)
	image_parser.add_argument("--height", type=int, required=True)
	image_parser.add_argument("--steps", type=int, required=True)
	args = parser.parse_args()
	evidence_path = args.evidence if args.command == "validate" else args.root.parent / (args.root.name + "-evidence.json")
	evidence = {"command": args.command, "job_id": os.environ.get("SLURM_JOB_ID"),
		"partition": os.environ.get("SLURM_JOB_PARTITION"), "passed": False}
	try:
		require(evidence["partition"] == "Debug" and evidence["job_id"], "checks require a Slurm Debug allocation")
		require(args.width > 0 and args.height > 0 and args.steps > 0, "invalid validation dimensions/steps")
		(validate if args.command == "validate" else images)(args, evidence)
		evidence["passed"] = True
	except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
		evidence["error"] = str(error)
		print(f"PFC_VALIDATION error={error}", file=sys.stderr)
	finally:
		evidence_path.parent.mkdir(parents=True, exist_ok=True)
		evidence_path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
	print(f"PFC_VALIDATION passed={str(evidence['passed']).lower()} evidence={evidence_path}")
	return 0 if evidence["passed"] else 1


if __name__ == "__main__":
	sys.exit(main())
