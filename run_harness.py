#!/usr/bin/env python3
"""CLI wrapper for the collision-time harness."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rh_engine.harness import main

if __name__ == "__main__":
    main()
