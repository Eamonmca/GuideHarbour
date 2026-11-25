"""Console entry point for the packaged Streamlit dashboard."""
from __future__ import annotations

from pathlib import Path
import sys

APP_FILENAME = "dashboard_demo.py"


def main() -> None:
    """Launch the current ProbeHarbour dashboard with Streamlit."""
    try:
        from streamlit.web import cli as streamlit_cli
    except ImportError as exc:
        raise SystemExit(
            "The dashboard dependencies are not installed. Run "
            "`pip install 'probeharbour[dashboard]'`."
        ) from exc
    # ``dashboard_demo.py`` is the current five-stage design workflow. Keep
    # the primary installed command pointed at it so a Git/wheel installation
    # opens the same interface as the repository development launcher.
    app = Path(__file__).with_name(APP_FILENAME)
    sys.argv = ["streamlit", "run", str(app), *sys.argv[1:]]
    raise SystemExit(streamlit_cli.main())


if __name__ == "__main__":
    main()
