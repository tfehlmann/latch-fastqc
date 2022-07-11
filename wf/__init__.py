"""
Process sequencing files with FASTQC.
"""

import os
import subprocess
from enum import Enum
from typing import List, Optional, Union

from latch import small_task, workflow
from latch.functions.messages import message
from latch.types import LatchAuthor, LatchDir, LatchFile, LatchMetadata, LatchParameter

from wf.helper import execute_cmd, latch2local


class Format(Enum):
    """
    Formats supported by FASTQC.
    """

    fastq = "fastq"  # pylint: disable=invalid-name
    bam = "bam"  # pylint: disable=invalid-name
    sam = "sam"  # pylint: disable=invalid-name
    bam_mapped = "bam_mapped"  # pylint: disable=invalid-name
    sam_mapped = "sam_mapped"  # pylint: disable=invalid-name


@small_task
def fastqc_task(  # pylint: disable=too-many-arguments
    input_files: List[Union[LatchFile, LatchDir]],
    output_folder: LatchDir,
    casava: bool,
    nano: bool,
    nofilter: bool,
    extract: bool,
    nogroup: bool,
    min_length: Optional[int],
    fmt: Optional[Format],
    contaminants: Optional[LatchFile],
    adapters: Optional[LatchFile],
    limits: Optional[LatchFile],
    kmers: int,
) -> LatchDir:
    """
    Run FASTQC command on input files
    """
    # pylint: disable=too-many-locals

    # get fastqc version
    cmd_version = ["fastqc", "--version"]
    result = subprocess.run(cmd_version, check=True, stdout=subprocess.PIPE, text=True)
    fastqc_version = result.stdout.split(" ")[-1]
    message(typ="info", data={"title": "FASTQC version", "body": fastqc_version})
    print("FASTQC version:", fastqc_version, flush=True)

    threads = str(len(os.sched_getaffinity(0)))
    fqc_output_folder = "/root/fastqc"
    os.makedirs(fqc_output_folder)
    cmd = [
        "fastqc",
        "--threads",
        threads,
        "--outdir",
        fqc_output_folder,
        "--casava" if casava else "",
        "--nano" if nano else "",
        "--nofilter" if nofilter else "",
        "--extract" if extract else "",
        "--nogroup" if nogroup else "",
        *(["--minlength", min_length] if min_length else []),
        *(["--format", fmt.value] if fmt else []),
        *(["--contaminants", contaminants] if contaminants else []),
        *(["--adapters", adapters] if adapters else []),
        *(["--limits", limits] if limits else []),
        *(["--kmers", kmers] if kmers else []),
    ]
    cmd.extend(latch2local(f) for f in input_files)
    # remove empty strings
    cmd = [e for e in cmd if e]
    message(typ="info", data={"title": "Running FASTQC", "body": f"Command: {' '.join(cmd)}"})
    execute_cmd(cmd, capture_stdout=True, capture_stderr=True)

    return LatchDir(
        path=fqc_output_folder,
        remote_path=output_folder.remote_path,  # type: ignore
    )


metadata = LatchMetadata(
    display_name="FASTQC",
    author=LatchAuthor(
        name="Tobias Fehlmann",
        email="tobias@astera.org",
        github="https://github.com/tfehlmann",
    ),
    repository="https://github.com/tfehlmann/latch-fastqc",
    license="MIT",
)

metadata.parameters["input_files"] = LatchParameter(
    display_name="Input files",
    description="Sequencing files to be processed. Must be FASTQ or BAM/SAM files.",
    hidden=False,
)

metadata.parameters["output_folder"] = LatchParameter(
    display_name="Output folder",
    description="Folder where FASTQC output files will be stored.",
    hidden=False,
    output=True,
)

metadata.parameters["casava"] = LatchParameter(
    display_name="casava",
    description="Files come from raw casava output. Files in the same sample "
    "group (differing only by the group number) will be analysed "
    "as a set rather than individually. Sequences with the filter "
    "flag set in the header will be excluded from the analysis. "
    "Files must have the same names given to them by casava "
    "(including being gzipped and ending with .gz) otherwise they "
    "won't be grouped together correctly.",
    hidden=False,
)

metadata.parameters["nano"] = LatchParameter(
    display_name="nano",
    description="Files come from nanopore sequences and are in fast5 format. In "
    "this mode you can pass in directories to process and the program "
    "will take in all fast5 files within those directories and produce "
    "a single output file from the sequences found in all files.",
    hidden=False,
)

metadata.parameters["nofilter"] = LatchParameter(
    display_name="nofilter",
    description="If running with --casava then don't remove read flagged by "
    "casava as poor quality when performing the QC analysis.",
    hidden=False,
)

metadata.parameters["extract"] = LatchParameter(
    display_name="extract",
    description="If set then the zipped output file will be uncompressed in "
    "the same directory after it has been created.",
    hidden=False,
)

metadata.parameters["nogroup"] = LatchParameter(
    display_name="nogroup",
    description="Disable grouping of bases for reads >50bp. All reports will "
    "show data for every base in the read.  WARNING: Using this "
    "option will cause fastqc to crash and burn if you use it on "
    "really long reads, and your plots may end up a ridiculous size. "
    "You have been warned!",
    hidden=False,
)

metadata.parameters["min_length"] = LatchParameter(
    display_name="min_length",
    description="Sets an artificial lower limit on the length of the sequence "
    "to be shown in the report. As long as you set this to a value "
    "greater or equal to your longest read length then this will be "
    "the sequence length used to create your read groups. This can "
    "be useful for making directly comaparable statistics from "
    "datasets with somewhat variable read lengths.",
    hidden=False,
)

metadata.parameters["format"] = LatchParameter(
    display_name="format",
    description="Bypasses the normal sequence file format detection and "
    "forces the program to use the specified format. Valid "
    "formats are bam,sam,bam_mapped,sam_mapped and fastq.",
    hidden=False,
)

metadata.parameters["contaminants"] = LatchParameter(
    display_name="contaminants",
    description="Specifies a non-default file which contains the list of "
    "contaminants to screen overrepresented sequences against. "
    "The file must contain sets of named contaminants in the "
    "form name[tab]sequence. Lines prefixed with a hash will be ignored.",
    hidden=False,
)

metadata.parameters["adapters"] = LatchParameter(
    display_name="adapters",
    description="Specifies a non-default file which contains the list of "
    "adapter sequences which will be explicity searched against "
    "the library. The file must contain sets of named adapters "
    "in the form name[tab]sequence. Lines prefixed with a hash will be ignored.",
    hidden=False,
)

metadata.parameters["limits"] = LatchParameter(
    display_name="limits",
    description="Specifies a non-default file which contains a set of criteria "
    "which will be used to determine the warn/error limits for the "
    "various modules. This file can also be used to selectively "
    "remove some modules from the output all together. The format "
    "needs to mirror the default limits.txt file found in the Configuration folder.",
    hidden=False,
)

metadata.parameters["kmers"] = LatchParameter(
    display_name="kmers",
    description="Specifies the length of Kmer to look for in the Kmer content "
    "module. Specified Kmer length must be between 2 and 10.",
    hidden=False,
)


@workflow(metadata=metadata)
def fastqc(  # pylint: disable=too-many-arguments
    input_files: List[Union[LatchFile, LatchDir]],
    output_folder: LatchDir,
    casava: bool = False,
    nano: bool = False,
    nofilter: bool = False,
    extract: bool = False,
    nogroup: bool = False,
    min_length: Optional[int] = None,
    format: Optional[Format] = None,  # pylint: disable=redefined-builtin
    contaminants: Optional[LatchFile] = None,
    adapters: Optional[LatchFile] = None,
    limits: Optional[LatchFile] = None,
    kmers: int = 7,
) -> LatchDir:
    """
    Run FASTQC command on input file(s)
    """

    return fastqc_task(
        input_files=input_files,
        output_folder=output_folder,
        casava=casava,
        nano=nano,
        nofilter=nofilter,
        extract=extract,
        nogroup=nogroup,
        min_length=min_length,
        fmt=format,
        contaminants=contaminants,
        adapters=adapters,
        limits=limits,
        kmers=kmers,
    )


if __name__ == "__main__":
    pass
