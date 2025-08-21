import subprocess


def quality_check():
    subprocess.run(
        "isort nnunet_contact_seg/*.py -c && black nnunet_contact_seg --check && snakefmt nnunet_contact_seg --check",
        shell=True,
        check=True,
    )


def quality_fix():
    subprocess.run(
        "isort nnunet_contact_seg/*.py && black nnunet_contact_seg && snakefmt nnunet_contact_seg",
        shell=True,
        check=True,
    )
