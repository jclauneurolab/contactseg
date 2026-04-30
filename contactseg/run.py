#!/usr/bin/env python3
import os
from pathlib import Path

from snakebids import bidsapp, plugins

# Import your utils module
try:
    from contactseg.workflow.lib import (
        utils_tf as utils_tf,
    )
except ImportError:
    from workflow.lib import utils_tf as utils_tf

if "__file__" not in globals():
    __file__ = "../contactseg/run.py"

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="contactseg"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)

# Set the conda prefix directory
conda_prefix = Path(utils_tf.get_cache_dir()) / "conda"

# Set templateflow directory to your OS cache via utils BEFORE the app runs
os.environ["TEMPLATEFLOW_HOME"] = str(Path(utils_tf.get_cache_dir()) / "templateflow")


# Set the environment variable SNAKEMAKE_CONDA_PREFIX if not already set
if "SNAKEMAKE_CONDA_PREFIX" not in os.environ:
    os.environ["SNAKEMAKE_CONDA_PREFIX"] = str(conda_prefix)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return app.build_parser().parser

if __name__ == "__main__":
    app.run()