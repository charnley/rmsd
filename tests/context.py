from pathlib import Path
import subprocess

RESOURCE_PATH = Path("tests/resources")

def call_main(args: list[str]) -> list[str]:

    root_path = Path("./")
    filename = root_path / "rmsd/calculate_rmsd.py"

    cmd = ["python", f"{filename}", *args]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()

    if stderr is not None:
        print(stderr.decode())

    return stdout.decode().strip().split("\n")
