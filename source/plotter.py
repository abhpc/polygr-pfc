#!/usr/bin/env python3
"""Export headerless PFC density fields as grayscale PNG images."""

import argparse
import os
from pathlib import Path
import re
import sys

import numpy as np
from PIL import Image


def safe_prefix(prefix):
	if not prefix:
		raise ValueError("output prefix must not be empty")
	return bytes(value if (65 <= value <= 90 or 97 <= value <= 122 or 48 <= value <= 57 or value in (45, 95))
		else 95 for value in os.fsencode(prefix)).decode("ascii")


def frame_step(path):
	match = re.fullmatch(r".+-t-([0-9]+)\.dat", Path(path).name)
	if match is None:
		match = re.fullmatch(r".+-([0-9]+)\.dat", Path(path).name)
	return int(match.group(1)) if match else None


def output_path(path, expected_step=None):
	path = Path(path)
	match = re.fullmatch(r"(.+)-([0-9]+)\.png", path.name)
	if match is None:
		raise ValueError(f"{path}: output must be <prefix>-<step>.png")
	step = int(match.group(2))
	if expected_step is not None and step != expected_step:
		raise ValueError(f"{path}: output step {step} differs from input step {expected_step}")
	return path.with_name(f"{safe_prefix(match.group(1))}-{step}.png")


def resolve_input(path):
	path = Path(path)
	try:
		path.stat()
		return path
	except FileNotFoundError:
		match = re.fullmatch(r"(.+)-t-([0-9]+)\.dat", path.name)
		if match is None:
			raise
		mapped = path.with_name(f"{safe_prefix(match.group(1))}-{int(match.group(2))}.dat")
		mapped.stat()
		print(f"PFC_INPUT requested={path} actual={mapped}")
		return mapped


def read_density(path, width, height):
	if width <= 0 or height <= 0:
		raise ValueError("field dimensions must be positive")
	data = np.loadtxt(path, dtype=np.float64, comments=None, ndmin=2)
	if data.shape != (width * height, 1):
		raise ValueError(f"{path}: expected {width * height} rows of one real value, found shape {data.shape}")
	if not np.isfinite(data).all():
		raise ValueError(f"{path}: density contains NaN or infinity")
	return data[:, 0]


def write_image(path, density, width, height):
	low, high = float(density.min()), float(density.max())
	if high == low:
		gray = np.zeros((height, width), dtype=np.uint8)
	else:
		span = high - low
		if np.isfinite(span):
			unit = (density - low) / span
		else:
			unit = (density * .5 - low * .5) / (high * .5 - low * .5)
		gray = np.flipud(np.floor(255 * unit + .5).astype(np.uint8).reshape(height, width))
	created = False
	try:
		with path.open("xb") as stream:
			created = True
			Image.fromarray(gray, mode="L").save(stream, format="PNG")
	except Exception:
		if created:
			path.unlink(missing_ok=True)
		raise


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input", metavar="INPUT")
	parser.add_argument("output", metavar="OUTPUT")
	parser.add_argument("width", metavar="W", type=int)
	parser.add_argument("height", metavar="H", type=int)
	parser.add_argument("interval", metavar="START INCREMENT END", nargs="*", type=int)
	args = parser.parse_args(argv)
	try:
		if args.width <= 0 or args.height <= 0:
			raise ValueError("W and H must be positive integers")
		if not args.input.endswith(".dat") or not args.output.endswith(".png"):
			raise ValueError("INPUT and OUTPUT must explicitly end in .dat and .png")
		batch = bool(args.interval)
		if batch:
			if len(args.interval) != 3 or args.input.count("#") != 1 or args.output.count("#") != 1:
				raise ValueError("batch mode needs START INCREMENT END and exactly one # in each template")
			start, increment, end = args.interval
			if start < 0 or increment <= 0 or end < start:
				raise ValueError("batch interval requires 0 <= START <= END and INCREMENT > 0")
			steps = range(start, end + 1, increment)
		else:
			if "#" in args.input or "#" in args.output:
				raise ValueError("single-frame paths must not contain #")
			steps = [None]
		targets = set()
		for step in steps:
			requested = Path(args.input.replace("#", str(step)) if batch else args.input)
			if batch and frame_step(requested) not in (None, step):
				raise ValueError(f"{requested}: template step differs from requested frame")
			output = args.output.replace("#", str(step)) if batch else args.output
			target = output_path(output, step if batch else frame_step(requested))
			identity = os.path.abspath(target)
			if identity in targets:
				raise ValueError(f"{target}: normalized output collision")
			if os.path.lexists(target):
				raise FileExistsError(f"{target}: refusing to overwrite existing output")
			targets.add(identity)
		generated, missing = 0, 0
		for step in steps:
			requested = Path(args.input.replace("#", str(step)) if batch else args.input)
			target = output_path(args.output.replace("#", str(step)) if batch else args.output,
				step if batch else frame_step(requested))
			try:
				source = resolve_input(requested)
				density = read_density(source, args.width, args.height)
			except FileNotFoundError:
				if not batch:
					raise
				print(f"Warning: {requested}: missing input frame", file=sys.stderr)
				missing += 1
				continue
			write_image(target, density, args.width, args.height)
			generated += 1
			print(f"PFC_PLOT_FRAME input={source} output={target}")
		if generated == 0:
			raise ValueError("no input frames were found")
		print(f"PFC_PLOT generated={generated} missing={missing}")
		return 0
	except (OSError, ValueError) as error:
		print(f"Error: {error}", file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
