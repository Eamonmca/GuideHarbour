"""Launch the guide designer packaged with the GuideHarbour fork."""
from pathlib import Path
import subprocess
import sys


def main():
    app = Path(__file__).with_name("guide_dashboard.py")
    raise SystemExit(subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), *sys.argv[1:]]))


if __name__ == "__main__":
    main()
