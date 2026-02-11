"""Generate all Jupyter notebooks for Chapters 10-11: Dynamic Game.

This master script generates and optionally executes all 5 notebooks:
  1. main_ch10_equilibrium.ipynb  - MPE computation and pseudo-data generation
  2. main_ch11_AM.ipynb           - Aguirregabiria-Mira (2007) estimation
  3. main_ch11_PSD.ipynb          - Pesendorfer-Schmidt-Dengler estimation
  4. main_ch11_forward_BBL.ipynb  - Forward Simulation P-SD + BBL inequality
  5. main_ch11_policy_sim.ipynb   - Counterfactual policy simulation

Usage:
  python generate_notebook_ch10_ch11.py           # Generate all notebooks
  python generate_notebook_ch10_ch11.py --execute  # Generate and execute all
  python generate_notebook_ch10_ch11.py --execute 1 5  # Generate all, execute #1 and #5 only

Execution order matters: notebooks 2-4 depend on Matlab data (data_from_matlab/),
notebook 5 is self-contained.
Estimated total execution time: ~2 hours (most time spent on BBL bootstrap in #4).
"""
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS = [
    ("generate_notebook_ch10_equilibrium.py", "main_ch10_equilibrium.ipynb"),
    ("generate_notebook_ch11_AM.py", "main_ch11_AM.ipynb"),
    ("generate_notebook_ch11_PSD.py", "main_ch11_PSD.ipynb"),
    ("generate_notebook_ch11_forward_BBL.py", "main_ch11_forward_BBL.ipynb"),
    ("generate_notebook_ch11_policy_sim.py", "main_ch11_policy_sim.ipynb"),
]

DESCRIPTIONS = [
    "Ch10: MPE computation and pseudo-data generation",
    "Ch11: Aguirregabiria-Mira (2007) estimation",
    "Ch11: Pesendorfer-Schmidt-Dengler estimation",
    "Ch11: Forward Simulation P-SD + BBL inequality estimator",
    "Ch11: Counterfactual policy simulation",
]


def generate_all():
    """Generate all notebooks by running each sub-generator."""
    script_dir = Path(__file__).parent
    print("=" * 70)
    print("Generating all notebooks for Chapters 10-11: Dynamic Game")
    print("=" * 70)

    for i, (script, notebook) in enumerate(SCRIPTS):
        print(f"\n[{i+1}/5] {DESCRIPTIONS[i]}")
        print(f"  Running: {script} -> {notebook}")
        result = subprocess.run(
            [sys.executable, str(script_dir / script)],
            cwd=str(script_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr}")
            return False
        print(f"  {result.stdout.strip()}")

    print("\n" + "=" * 70)
    print("All notebooks generated successfully.")
    print("=" * 70)
    return True


def execute_notebook(notebook_path):
    """Execute a notebook using jupyter nbconvert."""
    print(f"  Executing: {notebook_path.name} ...")
    start = time.time()
    result = subprocess.run(
        [
            sys.executable, "-m", "jupyter", "nbconvert",
            "--to", "notebook",
            "--execute",
            "--ExecutePreprocessor.timeout=7200",
            str(notebook_path),
        ],
        capture_output=True,
        text=True,
    )
    elapsed = time.time() - start
    if result.returncode != 0:
        print(f"  FAILED ({elapsed:.0f}s): {result.stderr[:200]}")
        return False
    print(f"  Done ({elapsed:.0f}s)")
    return True


def main():
    args = sys.argv[1:]
    do_execute = "--execute" in args
    execute_indices = set()

    if do_execute:
        args.remove("--execute")
        if args:
            execute_indices = {int(x) for x in args}
        else:
            execute_indices = {1, 2, 3, 4, 5}

    # Step 1: Generate all notebooks
    if not generate_all():
        sys.exit(1)

    # Step 2: Execute if requested
    if do_execute:
        script_dir = Path(__file__).parent
        print("\n" + "=" * 70)
        print("Executing notebooks")
        print("=" * 70)

        for i, (_, notebook) in enumerate(SCRIPTS):
            if (i + 1) in execute_indices:
                nb_path = script_dir / notebook
                print(f"\n[{i+1}/5] {DESCRIPTIONS[i]}")
                execute_notebook(nb_path)
            else:
                print(f"\n[{i+1}/5] Skipped: {notebook}")

        print("\n" + "=" * 70)
        print("Execution complete.")
        print("=" * 70)

    # Summary
    script_dir = Path(__file__).parent
    print("\nGenerated notebooks:")
    for i, (_, notebook) in enumerate(SCRIPTS):
        nb_path = script_dir / notebook
        status = "exists" if nb_path.exists() else "MISSING"
        print(f"  {i+1}. {notebook} [{status}]")

    print(f"\nOutput directory: {(script_dir / '..' / 'output').resolve()}")


if __name__ == "__main__":
    main()
