#!/usr/bin/env python

import argparse
import bz2
import gzip
import lzma
import math
import os
import re
import subprocess
import sys
import tempfile
from enum import Enum, auto
from concurrent.futures import ThreadPoolExecutor
from shutil import which


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    noncompressed = auto()


def is_compressed(file):
    with open(file, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.noncompressed


def count_sequences(file):
    n = 0
    compression = is_compressed(file)
    if compression == Compression.gzip:
        fasta = gzip.open(file, "rt")
    elif compression == Compression.bzip2:
        fasta = bz2.open(file, "rt")
    elif compression == Compression.xz:
        fasta = lzma.open(file, "rt")
    else:
        fasta = open(file, "r")
    for line in fasta:
        if line.startswith(">"):
            n += 1
    fasta.close()
    return n


def run_prodigal(opts, workDir, currentId, chunkFile):
    cmd = ["prodigal-gv", "-q", "-p", "meta", "-i", chunkFile.name]

    if opts.proteins:
        cmd.append("-a")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".faa")

    if opts.closed:
        cmd.append("-c")

    if opts.nucl:
        cmd.append("-d")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".fna")

    if opts.format:
        cmd.append("-f")
        cmd.append(opts.format)

    if opts.mask:
        cmd.append("-m")

    if opts.scorefile:
        cmd.append("-s")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".score")

    cmd.append("-o")
    cmd.append(workDir.name + "/chunk" + str(currentId) + ".out")

    subprocess.run(cmd, shell=False, check=True)


def append_fasta_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line.startswith(">"):
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = (
                        match.group(1)
                        + str(startNum)
                        + "_"
                        + match.group(3)
                        + match.group(4)
                    )
                trgt.write(line.strip() + "\n")


def append_gff_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line.startswith("##gff-version"):
                    continue
                elif line[0] != "#" and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = (
                        match.group(1)
                        + str(startNum)
                        + "_"
                        + match.group(3)
                        + match.group(4)
                    )
                trgt.write(line.strip() + "\n")


def append_gbk_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line[0] == " " and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = (
                        match.group(1)
                        + str(startNum)
                        + "_"
                        + match.group(3)
                        + match.group(4)
                    )
                trgt.write(line.strip() + "\n")


def append_raw_file(file, targetFile):
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            trgt.write(input.read())


def print_gff_file(file, startNum):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            if line.startswith("##gff-version"):
                continue
            elif line[0] != "#" and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = (
                    match.group(1)
                    + str(startNum)
                    + "_"
                    + match.group(3)
                    + match.group(4)
                )
            print(line.strip() + "\n")


def print_gbk_file(file, startNum):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            if line[0] == " " and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = (
                    match.group(1)
                    + str(startNum)
                    + "_"
                    + match.group(3)
                    + match.group(4)
                )
            print(line.strip() + "\n")


def print_raw_file(file):
    with open(file, "r") as input:
        trgt.write(input.read())

def parse_cli():
    argp = argparse.ArgumentParser(description="Parallel prodigal-gv gene prediction")
    argp.add_argument(
        "-i",
        "--input",
        type=str,
        help="Specify FASTA/Genbank input file (default reads from stdin).",
    )
    argp.add_argument(
        "-d",
        "--nucl",
        type=str,
        help="Write nucleotide sequences of genes to the selected file.",
    )
    argp.add_argument(
        "-a",
        "--proteins",
        type=str,
        help="Write protein translations to the selected file.",
    )
    argp.add_argument(
        "-o",
        "--output",
        type=str,
        help="Specify output file (default writes to stdout).",
    )
    argp.add_argument(
        "-s",
        "--scorefile",
        type=str,
        help="Write all potential genes (with scores) to the selected file.",
    )
    argp.add_argument(
        "-f",
        "--format",
        type=str,
        default="gbk",
        help="Select output format (gbk, gff, or sco). Default is gbk.",
    )
    argp.add_argument(
        "-c",
        "--closed",
        action="store_true",
        help="Closed ends. Do not allow genes to run off edges.",
    )
    argp.add_argument(
        "-m",
        "--mask",
        action="store_true",
        help="Treat runs of N as masked sequence; don't build genes across them.",
    )
    argp.add_argument(
        "-t",
        "--tasks",
        type=int,
        default=os.cpu_count(),
        help="number of prodigal-gv processes to start in parallel.",
    )
    argp.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Run quietly (will not print output to stout).",
    )
    opts = argp.parse_args()
    return argp, opts

def main():
    argp, opts = parse_cli()

    if len(sys.argv) == 1:
        argp.print_help()
        sys.exit(0)

    if opts.tasks is not None:
        opts.tasks = max(opts.tasks, 1)
        tasks = opts.tasks

    if which("prodigal-gv") is None:
        raise ValueError("prodigal-gv not found in the PATH.")

    # Delete opput files if they already exist. Otherwise the new results will
    # be appended to them
    if opts.output and os.path.isfile(opts.output):
        os.remove(opts.output)
    if opts.nucl and os.path.isfile(opts.nucl):
        os.remove(opts.nucl)
    if opts.proteins and os.path.isfile(opts.proteins):
        os.remove(opts.proteins)
    if opts.scorefile and os.path.isfile(opts.scorefile):
        os.remove(opts.scorefile)

    seqCnt = 0
    currentChunk = 1

    workDir = tempfile.TemporaryDirectory()
    executor = ThreadPoolExecutor(max_workers=tasks)
    currentFile = open(workDir.name + "/chunk" + str(currentChunk), "w")

    queryFile = None
    if not opts.input:
        queryFile = "/dev/fd/0"
        if sys.stdin.isatty():
            print("Cannot read sequences from STDIN.")
            exit(1)
    else:
        queryFile = opts.input

    n_sequences = count_sequences(queryFile)
    seqsPerChunk = math.ceil(n_sequences / tasks)

    # Determine if the file is compressed and open it.
    compression = is_compressed(queryFile)
    if compression == Compression.gzip:
        fasta = gzip.open(queryFile, "rt")
    elif compression == Compression.bzip2:
        fasta = bz2.open(queryFile, "rt")
    elif compression == Compression.xz:
        fasta = lzma.open(queryFile, "rt")
    else:
        fasta = open(queryFile, "r")

    for line in fasta:
        if line[0] == ">" and seqCnt == seqsPerChunk:
            currentFile.close()
            executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)
            currentFile = None
            seqCnt = 0
            currentChunk += 1
        if currentFile is None:
            currentFile = open(workDir.name + "/chunk" + str(currentChunk), "w")
        currentFile.write(line)
        if line[0] == ">":
            seqCnt += 1

    fasta.close()

    if seqCnt > 0:
        currentFile.close()
        executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)

    # await completion of tasks
    executor.shutdown(wait=True)

    # collect output
    proteinFile = opts.proteins
    nuclFile = opts.nucl
    outFile = opts.output
    scoreFile = opts.scorefile

    protIdStart = 0
    nuclIdStart = 0
    gffIdStart = 0
    gbkIdStart = 0
    for cur in range(1, currentChunk + 1):
        if proteinFile:
            append_fasta_file(
                workDir.name + "/chunk" + str(cur) + ".faa", protIdStart, proteinFile
            )
        if nuclFile:
            append_fasta_file(
                workDir.name + "/chunk" + str(cur) + ".fna", nuclIdStart, nuclFile
            )
        if scoreFile:
            append_raw_file(workDir.name + "/chunk" + str(cur) + ".score", scoreFile)

        if outFile:
            if opts.format == "gbk":
                append_gbk_file(
                    workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart, outFile
                )
            elif opts.format == "gff":
                append_gff_file(
                    workDir.name + "/chunk" + str(cur) + ".out", gffIdStart, outFile
                )
            else:
                append_raw_file(workDir.name + "/chunk" + str(cur) + ".out", outFile)

        elif not opts.quiet:
            if opts.format == "gbk":
                print_gbk_file(
                    workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart
                )
            elif opts.format == "gff":
                print_gff_file(
                    workDir.name + "/chunk" + str(cur) + ".out", gffIdStart
                )
            else:
                print_raw_file(workDir.name + "/chunk" + str(cur) + ".out")

        protIdStart += seqsPerChunk
        nuclIdStart += seqsPerChunk
        gffIdStart += seqsPerChunk
        gbkIdStart += seqsPerChunk


if __name__ == "__main__":
    main()
