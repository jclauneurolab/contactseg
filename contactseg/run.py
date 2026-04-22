#!/usr/bin/env python3
import os
from pathlib import Path

from snakebids import bidsapp, plugins

# Import your utils module
try:
    from contactseg.workflow.lib import utils_tf as utils_tf
except ImportError:
    from workflow.lib import utils_tf as utils_tf

if "__file__" not in globals():
    __file__ = "../contactseg/run.py"

# Set templateflow directory to your OS cache via utils BEFORE the app runs
os.environ["TEMPLATEFLOW_HOME"] = str(Path(utils_tf.get_download_dir()) / "templateflow")

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="contactseg"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)

def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return app.build_parser().parser

if __name__ == "__main__":
    app.run()