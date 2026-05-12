"""Package imports and CLI loads."""

from typer.testing import CliRunner

import portein
from portein.__main__ import app


def test_public_api_present():
    for name in portein.__all__:
        assert hasattr(portein, name), f"portein.{name} missing"


def test_cli_lists_all_subcommands():
    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    for sub in ("rotate", "secondary", "pymol", "illustrate"):
        assert sub in result.output
