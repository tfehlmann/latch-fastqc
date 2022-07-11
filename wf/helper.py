import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Union

from latch.functions.messages import message
from latch.types import LatchDir, LatchFile


def execute_cmd(cmd, capture_stdout=True, capture_stderr=True):
    """
    Execute a command in a shell and print captured output.
    """
    env = dict(os.environ)
    env["PATH"] = f"/root/micromamba/bin:{env['PATH']}"
    cmd_str = shlex.join(cmd)
    # unquote stream redirects
    cmd_str = (
        cmd_str.replace(" '>' ", " > ")
        .replace(" '<' ", " < ")
        .replace(" '2>&1' ", " 2>&1 ")
        .replace(" '|' ", " | ")
    )
    print(f"Execute cmd: {cmd_str}", flush=True)
    stdout = subprocess.PIPE if capture_stdout else None
    stderr = subprocess.PIPE if capture_stderr else None
    if capture_stdout and capture_stderr:
        stderr = subprocess.STDOUT

    with subprocess.Popen(
        cmd_str,
        stdout=stdout,
        stderr=stderr,
        env=env,
        executable="/bin/bash",
        shell=True,
    ) as process:
        if capture_stdout:
            for line in process.stdout:  # type: ignore
                print(line.decode("utf-8"))
        elif capture_stderr:
            for line in process.stderr:  # type: ignore
                print(line.decode("utf-8"), file=sys.stderr)
                message(
                    typ="warning",
                    data={
                        "title": "Warning/Error while running FASTQC",
                        "body": line.decode("utf-8"),
                    },
                )
        status = process.wait()
        if status != 0:
            print(f"Failure! Command exited with status {status}.", file=sys.stderr)
            message(
                typ="error",
                data={
                    "title": "Error while running FASTQC",
                    "body": f"Failure! Command exited with status {status}.",
                },
            )
            raise subprocess.CalledProcessError(status, cmd=cmd)


def latch2local(latchdir_or_file: Union[LatchFile, LatchDir]):
    """
    Convert a latch file or directory to a local file or directory.
    """
    return str(Path(latchdir_or_file).resolve())
