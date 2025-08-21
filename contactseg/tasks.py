import subprocess


def quality_check():
    subprocess.run(
        "isort contactseg/*.py -c && black contactseg --check && snakefmt contactseg --check",
        shell=True,
        check=True,
    )


def quality_fix():
    subprocess.run(
        "isort contactseg/*.py && black contactseg && snakefmt contactseg",
        shell=True,
        check=True,
    )
