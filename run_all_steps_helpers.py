"""Helpers shared by run_all_steps_to_generate_the_hg38_catalog.py and run_all_steps_to_generate_the_CHM13_catalog.py.

run() and download_unless_local_copy_is_up_to_date() take the pipeline's parsed command-line args, which must include
dry_run, only_step, start_with_step and end_with_step. Each pipeline script binds its args with functools.partial,
so its call sites can stay run(command, step_number=N).
"""

import gzip
import json
import os
import re
import requests
import subprocess
import sys
import tempfile


def run(command, step_number=None, *, args):
    """Runs a shell command, unless --dry-run or the --only-step, --start-with-step or --end-with-step range skips it.

    Args:
        command (str): shell command
        step_number (int): pipeline step that this command belongs to. Commands without a step number always run
            (except in a dry run).
        args (argparse.Namespace): the pipeline's parsed command-line args
    """
    command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
    if args.dry_run:
        if step_number is not None:
            print(f"STEP #{step_number}: {command}")
        else:
            print(command)

    if not args.dry_run or command.startswith("mkdir"):
        if step_number is None or (
            not (args.only_step is not None and step_number != args.only_step) and
            not (args.start_with_step is not None and step_number < args.start_with_step) and
            not (args.end_with_step is not None and step_number > args.end_with_step)
        ):
            if step_number is not None:
                print(f"\nSTEP #{step_number}: {command}")
            else:
                print(f"\n--- $ {command}")
            subprocess.run(command, shell=True, check=True)

def check_that_merge_loci_ignores_discarded_loci():
    """Raises RuntimeError unless the installed str_analysis.merge_loci keeps a locus that only matches a locus it
    already discarded.

    A merge_loci from before that str-analysis fix still compares each new locus with the loci it discarded, and drops
    the new locus if one of them matches it, although nothing in the output describes it. The pipelines' merges then
    silently lose loci (3,185 of them in the hg38 v2.1 catalog of 10/19/2026). This merges three (A)* loci with the
    pipelines' keep-first and 0.66 Jaccard similarity settings: the second is redundant with the first (Jaccard
    similarity 18/24), and the third only with the second (24/34, but 18/34 with the first), so the first and the third
    should be kept.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        catalog_paths = []
        for i, region in enumerate(["chr1:1000-1018", "chr1:1000-1024", "chr1:1000-1034"]):
            catalog_paths.append(os.path.join(temp_dir, f"catalog{i}.json"))
            with open(catalog_paths[-1], "wt") as f:
                json.dump([{"LocusId": f"locus{i}", "LocusStructure": "(A)*", "ReferenceRegion": region,
                            "VariantType": "Repeat"}], f)
        output_prefix = os.path.join(temp_dir, "merged")
        subprocess.run([
            sys.executable, "-m", "str_analysis.merge_loci", "--output-format", "JSON",
            "--overlapping-loci-action", "keep-first", "--overlap-fraction", "0.66", "--min-jaccard-similarity", "0.66",
            "--output-prefix", output_prefix, *catalog_paths], check=True, capture_output=True)
        with gzip.open(f"{output_prefix}.json.gz", "rt") as f:
            locus_ids = [record["LocusId"] for record in json.load(f)]

    if locus_ids != ["locus0", "locus2"]:
        raise RuntimeError(f"The installed str_analysis.merge_loci kept {locus_ids} rather than ['locus0', 'locus2'], "
                           f"so it lets loci it discarded make later loci redundant. Upgrade it to a version that "
                           f"includes str-analysis commit 383709c.")

def chdir(d):
    print(f"cd {d}")
    os.chdir(d)

def is_bgzipped(path):
    """Returns True if path is bgzipped, and so can be tabix-indexed, rather than plain gzipped.

    A bgzipped file starts with a gzip header that sets the FEXTRA flag and carries a "BC"
    extra subfield. Checking the bytes means this works for a file that a previous run
    already compressed, where the branch that compressed it is no longer known.
    """
    with open(path, "rb") as f:
        header = f.read(14)
    return header[:4] == b"\x1f\x8b\x08\x04" and header[12:14] == b"BC"

# The LongTR, HipSTR and GangSTR converters in str-analysis write start_0based + 1, so those three
# formats are 1-based rather than true BED. Indexing them with "-p bed", which means "-s1 -b2 -e3 -0",
# shifts every record one base to the right and makes tabix miss it at the left edge of a query.
ONE_BASED_BED_SUFFIXES = (".LongTR.bed.gz", ".HipSTR.bed.gz", ".GangSTR.bed.gz")

def tabix_coordinate_flags(path):
    """Returns the tabix flags describing path's coordinate convention."""
    return "-s1 -b2 -e3" if path.endswith(ONE_BASED_BED_SUFFIXES) else "-p bed"

def get_remote_file_etag(url):
    """Returns the ETag that the server reports for url. It changes whenever the file's contents change."""
    if url.startswith("gs://"):
        return subprocess.run(["gcloud", "storage", "objects", "describe", url, "--format=value(etag)"],
                              capture_output=True, text=True, check=True).stdout.strip()

    response = requests.head(url, allow_redirects=True, timeout=60)
    response.raise_for_status()
    if not response.headers.get("ETag"):
        raise ValueError(f"{url} has no ETag, so there's no way to tell whether a local copy of it is up to date")
    return response.headers["ETag"]

def download_unless_local_copy_is_up_to_date(url, local_dir, *, args):
    """Downloads url into local_dir, unless the copy already there has the file's current ETag.

    The ETag of each download is saved next to it in a .etag file, which is what later runs compare against.

    Args:
        url (str): gs:// path or URL of the file
        local_dir (str): directory to download the file into
        args (argparse.Namespace): the pipeline's parsed command-line args

    Returns:
        str: path of the local copy
    """
    local_path = os.path.join(local_dir, os.path.basename(url))
    if args.dry_run:
        print(f"Would download {url} to {local_path} unless the local copy is up to date")
        return local_path

    etag_path = f"{local_path}.etag"
    remote_etag = get_remote_file_etag(url)
    if os.path.isfile(local_path) and os.path.isfile(etag_path):
        with open(etag_path) as f:
            if f.read() == remote_etag:
                print(f"Local copy of {url} is up to date: {local_path}")
                return local_path

    if url.startswith("gs://"):
        run(f"gcloud storage cp {url} {local_path}.tmp", args=args)
    else:
        run(f"wget -q -O {local_path}.tmp {url}", args=args)
    run(f"mv {local_path}.tmp {local_path}", args=args)
    with open(etag_path, "w") as f:
        f.write(remote_etag)

    return local_path
