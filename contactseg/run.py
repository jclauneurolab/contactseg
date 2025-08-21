#!/usr/bin/env python3
import os
import sys
from pathlib import Path

from snakebids import bidsapp, plugins

if "__file__" not in globals():
    __file__ = "../contactseg/run.py"
import os

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
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
