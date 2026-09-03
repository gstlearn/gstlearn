#!/usr/bin/env python3

"""
===============================================================================
GitHub Workflow Execution Time Extractor (Demo Tests)
===============================================================================

Usage:
  python3 demo_times.py FIRST_RUN LAST_RUN [TEST_NAME]

Prerequisites:
  You must set the GITHUB_TOKEN environment variable before running the script:
    csh/tcsh : setenv GITHUB_TOKEN ghp_YourGitHubTokenHere
    bash/zsh : export GITHUB_TOKEN=ghp_YourGitHubTokenHere
  This can also be done by simply running the script and following the instructions:
    python3 get_token.py

Examples:
  python3 demo_times.py 3510 3510 Tuto_2D_ipynb
===============================================================================
"""

import os
import sys
import json
import re
import io
import zipfile
import urllib.request
import urllib.error
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

OWNER = "gstlearn"
REPO = "gstlearn"
WORKFLOW_NAME = "nonreg-demos-courses_ubuntu-latest"

GITHUB_TOKEN = os.getenv("GITHUB_TOKEN", None)


# ============================================================================
# HELP & USAGE
# ============================================================================


def print_help(error_msg=None):
    if error_msg:
        print(f"Error: {error_msg}\n", file=sys.stderr)

    help_text = f"""Usage: python3 {sys.argv[0]} FIRST_RUN LAST_RUN [TEST_NAME]

Fetch and display execution times for successful '{WORKFLOW_NAME}' workflow runs.

Arguments:
  FIRST_RUN     Positive integer representing the start run number.
  LAST_RUN      Positive integer representing the end run number.
  TEST_NAME     (Optional) Specific test name to parse from logs (e.g. Tuto_2D_ipynb).

Environment Variable:
  GITHUB_TOKEN  Required to avoid API rate limits and to download log archives.

Examples:
  python3 {sys.argv[0]} 3510 3510 Tuto_2D_ipynb
"""
    print(help_text, file=sys.stderr)


# ============================================================================
# GITHUB API & LOGS
# ============================================================================


class NoRedirectionHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def get_raw_response(url, auth=True):
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "gstlearn-demo-times",
    }
    if auth and GITHUB_TOKEN:
        headers["Authorization"] = f"Bearer {GITHUB_TOKEN}"

    req = urllib.request.Request(url, headers=headers)
    try:
        return urllib.request.urlopen(req)
    except urllib.error.HTTPError as e:
        if e.code == 403:
            print(
                f"\n[API Error 403] GitHub rate limit reached or access forbidden!",
                file=sys.stderr,
            )
        else:
            print(f"\n[API Error {e.code}] {e.reason}", file=sys.stderr)
        raise e
    except Exception as e:
        print(f"\n[Network Error] {e}", file=sys.stderr)
        raise e


def get_json(url):
    with get_raw_response(url) as response:
        return json.load(response)


def get_runs(first_run):
    runs = []
    page = 1

    while True:
        url = (
            f"https://api.github.com/repos/{OWNER}/{REPO}"
            f"/actions/workflows/{WORKFLOW_NAME}.yml/runs"
            f"?per_page=100&page={page}"
        )

        try:
            data = get_json(url)
        except Exception:
            break

        page_runs = data.get("workflow_runs", [])
        if not page_runs:
            break

        runs.extend(page_runs)

        oldest_in_page = page_runs[-1]["run_number"]
        if oldest_in_page <= first_run:
            break

        page += 1

    return runs


def identify_platform(filepath, content):
    """Identify target platform/job (python, r) from path or log content."""
    path_lower = filepath.lower()
    content_lower = content.lower()

    if "python" in path_lower or "call-python" in path_lower:
        return "python"
    if "/r/" in path_lower or "call-r" in path_lower:
        return "r"

    if "execute python" in content_lower or "pytest" in content_lower:
        return "python"
    if "execute r" in content_lower or "rscript" in content_lower:
        return "r"

    return None


def get_test_times_from_logs(run_id, test_name):
    url = f"https://api.github.com/repos/{OWNER}/{REPO}/actions/runs/{run_id}/logs"
    results = {}

    opener = urllib.request.build_opener(NoRedirectionHandler())
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "gstlearn-demo-times",
    }
    if GITHUB_TOKEN:
        headers["Authorization"] = f"Bearer {GITHUB_TOKEN}"

    req = urllib.request.Request(url, headers=headers)

    download_url = None
    try:
        with opener.open(req) as response:
            pass
    except urllib.error.HTTPError as e:
        if e.code in (301, 302, 307):
            download_url = e.headers.get("Location")
        else:
            return results
    except Exception:
        return results

    if not download_url:
        return results

    try:
        with get_raw_response(download_url, auth=False) as response:
            zip_bytes = response.read()
    except Exception:
        return results

    # Regex tolérante : accepte les points et les underscores dans le nom transmis
    normalized_test = re.escape(test_name.replace(".", "_"))
    pattern = re.compile(
        r"Test\s+#\d+:\s+" + normalized_test + r"\b[^\n\r]*?([\d\.]+)\s+sec",
        re.IGNORECASE,
    )

    try:
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
            for filename in zf.namelist():
                if not filename.endswith(".txt"):
                    continue

                with zf.open(filename) as log_file:
                    content = log_file.read().decode("utf-8", errors="ignore")

                    match = pattern.search(content)
                    if match:
                        key = identify_platform(filename, content)
                        if key and key not in results:
                            results[key] = float(match.group(1))
    except Exception:
        pass

    return results


# ============================================================================
# TIME UTILS
# ============================================================================


def elapsed_seconds(started_at, completed_at):
    if not started_at or not completed_at:
        return None
    t0 = datetime.fromisoformat(started_at.replace("Z", "+00:00"))
    t1 = datetime.fromisoformat(completed_at.replace("Z", "+00:00"))
    return (t1 - t0).total_seconds()


def format_elapsed(seconds):
    if seconds is None:
        return "-"
    seconds = int(round(seconds))
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60

    if hours:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    return f"{minutes}m {secs:02d}s"


def format_run_date(created_at):
    if not created_at:
        return "-"
    dt = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
    return dt.strftime("%Y-%m-%d")


# ============================================================================
# RUN JOBS
# ============================================================================


def get_build_jobs(run_id):
    url = f"https://api.github.com/repos/{OWNER}/{REPO}/actions/runs/{run_id}/jobs?per_page=100"
    try:
        data = get_json(url)
    except Exception:
        return {}

    result = {}
    for job in data.get("jobs", []):
        name = job.get("name", "").lower()
        key = None
        if "call-python" in name:
            key = "python"
        elif "call-r" in name:
            key = "r"

        if key:
            result[key] = {
                "elapsed": elapsed_seconds(
                    job.get("started_at"), job.get("completed_at")
                ),
            }
    return result


# ============================================================================
# CHART GENERATION
# ============================================================================


def generate_chart(history_data, test_name=None):
    """Generates 2 stacked subplots (Python and R) with independent Y-scales."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "\n[Warning] 'matplotlib' is not installed. PNG plot skipped.",
            file=sys.stderr,
        )
        print("Install it via: pip install matplotlib\n", file=sys.stderr)
        return

    if not history_data:
        return

    runs = [d["run"] for d in history_data]

    def to_min(val_sec):
        return val_sec / 60.0 if val_sec is not None else None

    environments = [
        ("Python", [to_min(d["python"]) for d in history_data], "#3572A5"),
        ("R", [to_min(d["r"]) for d in history_data], "#198CE7"),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    title_suffix = f" (Test: {test_name})" if test_name else " (Job Total)"
    fig.suptitle(
        f"gstlearn Demos - Execution Times{title_suffix}",
        fontsize=14,
        fontweight="bold",
    )

    for i, (name, vals, color) in enumerate(environments):
        ax = axes[i]
        ax.plot(runs, vals, marker="o", label=name, color=color, linewidth=2)
        ax.set_ylabel("Duration (min)", fontsize=9)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(loc="upper right", fontsize=9)

        # Lignes verticales à chaque changement de date
        last_date = None
        for item in history_data:
            current_date = item["date"]
            current_run = item["run"]
            if last_date is None or current_date != last_date:
                ax.axvline(
                    x=current_run, color="gray", linestyle=":", alpha=0.5, linewidth=1
                )
                if i == 0:  # Affiche la date seulement sur le premier subplot
                    ax.text(
                        current_run,
                        ax.get_ylim()[1],
                        f" {current_date}",
                        rotation=90,
                        verticalalignment="top",
                        horizontalalignment="left",
                        fontsize=7,
                        color="gray",
                    )
                last_date = current_date

    axes[-1].set_xlabel("Workflow Run Number", fontsize=11)
    plt.tight_layout()

    filename = f"demo_times_{test_name}.png" if test_name else "demo_times.png"
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Chart saved to {filename}")


# ============================================================================
# MAIN PROGRAM
# ============================================================================


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help"):
        print_help()
        sys.exit(0)

    if len(sys.argv) < 3:
        print_help("Missing required run range arguments.")
        sys.exit(1)

    try:
        first_run = int(sys.argv[1])
        last_run = int(sys.argv[2])
    except ValueError:
        print_help("FIRST_RUN and LAST_RUN must be positive integers.")
        sys.exit(1)

    if first_run <= 0 or last_run <= 0:
        print_help("Run numbers must be strictly greater than 0.")
        sys.exit(1)

    test_name = sys.argv[3] if len(sys.argv) >= 4 else None

    if first_run > last_run:
        first_run, last_run = last_run, first_run

    target_info = f" (Test: {test_name})" if test_name else " (Job Total)"
    print(
        f"\nSearching for workflows {first_run} to {last_run}{target_info}...",
        flush=True,
    )

    runs = get_runs(first_run)
    runs_in_order = [
        r
        for r in runs
        if r.get("run_number") and first_run <= r["run_number"] <= last_run
    ]
    runs_in_order.sort(key=lambda x: x["run_number"])

    if not runs_in_order:
        print("No matching runs found (or API rate limit exceeded).", flush=True)
        sys.exit(0)

    headers = ["Run", "Date", "Python", "R"]
    widths = [8, 13, 14, 14]

    print("\n" + "".join(f"{h:<{w}}" for h, w in zip(headers, widths)), flush=True)
    print("-" * sum(widths), flush=True)

    history_data = []

    for run in runs_in_order:
        if run.get("conclusion") != "success":
            continue

        internal_id = run["id"]
        run_date = format_run_date(run.get("created_at"))

        if test_name:
            test_times = get_test_times_from_logs(internal_id, test_name)
            python_sec = test_times.get("python")
            r_sec = test_times.get("r")

            python_val = f"{python_sec:.2f}s" if python_sec is not None else "-"
            r_val = f"{r_sec:.2f}s" if r_sec is not None else "-"
        else:
            jobs = get_build_jobs(internal_id)
            python_sec = jobs.get("python", {}).get("elapsed")
            r_sec = jobs.get("r", {}).get("elapsed")

            python_val = format_elapsed(python_sec)
            r_val = format_elapsed(r_sec)

        history_data.append(
            {
                "run": run["run_number"],
                "date": run_date,
                "python": python_sec,
                "r": r_sec,
            }
        )

        values = [
            str(run["run_number"]),
            run_date,
            python_val,
            r_val,
        ]

        print("".join(f"{v:<{w}}" for v, w in zip(values, widths)), flush=True)

    print()
    generate_chart(history_data, test_name)


if __name__ == "__main__":
    main()
