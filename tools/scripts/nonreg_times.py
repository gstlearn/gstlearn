#!/usr/bin/env python3

"""
===============================================================================
GitHub Workflow Execution Time Extractor (with optional test-level detail)
===============================================================================

Usage:
  python3 nonreg_times.py FIRST_RUN LAST_RUN [TEST_NAME]

Prerequisites:
  You must set the GITHUB_TOKEN environment variable before running the script:
    csh/tcsh : setenv GITHUB_TOKEN ghp_YourGitHubTokenHere
    bash/zsh : export GITHUB_TOKEN=ghp_YourGitHubTokenHere
  This can also be done by simply running the script and following the instructions:
    python3 get_token.py

Examples:
  python3 nonreg_times.py 3900 3936
  python3 nonreg_times.py 3936 3936 bench_Migrate
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
WORKFLOW_NAME = "nonreg-tests"

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
  TEST_NAME     (Optional) Specific test name to parse from logs (e.g. bench_Migrate).

Environment Variable:
  GITHUB_TOKEN  Required to avoid API rate limits and to download log archives.
                Set it in your shell prior to execution:
                  csh/tcsh : setenv GITHUB_TOKEN ghp_YourTokenHere
                  bash/zsh : export GITHUB_TOKEN=ghp_YourTokenHere

Examples:
  python3 {sys.argv[0]} 3900 3936
  python3 {sys.argv[0]} 3936 3936 bench_Migrate
"""
    print(help_text, file=sys.stderr)


# ============================================================================
# GITHUB API & LOGS
# ============================================================================


class NoRedirectionHandler(urllib.request.HTTPRedirectHandler):
    """Prevents urllib from automatically following 302 redirects."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def get_raw_response(url, auth=True):
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "gstlearn-nonreg-times",
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
            print(
                "Ensure GITHUB_TOKEN is properly exported in your environment (setenv GITHUB_TOKEN ghp_...).",
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
    """Identify target platform (msvc, ubuntu, rtools, macos) from path or content."""
    lower_path = filepath.lower()

    if "msvc" in lower_path:
        return "msvc"
    if "rtools" in lower_path:
        return "rtools"
    if "ubuntu" in lower_path:
        return "ubuntu"
    if "macos" in lower_path or "mac" in lower_path:
        return "macos"

    if "msvc" in content or "cl.exe" in content:
        return "msvc"
    if "rtools" in content or "mingw" in content:
        return "rtools"
    if "ubuntu" in content or "linux" in content:
        return "ubuntu"
    if "macos" in content or "darwin" in content or "apple" in content:
        return "macos"

    return None


def get_test_times_from_logs(run_id, test_name):
    """Downloads the run logs ZIP archive using the internal API run_id and extracts test durations in seconds."""
    url = f"https://api.github.com/repos/{OWNER}/{REPO}/actions/runs/{run_id}/logs"
    results = {}

    opener = urllib.request.build_opener(NoRedirectionHandler())
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "gstlearn-nonreg-times",
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

    pattern = re.compile(
        re.escape(test_name)
        + r"\b[^\n\r]*?(?:Passed|Completed|Took)[^\n\r]*?([\d\.]+)\s+sec",
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
        name = job.get("name", "")
        if not name.endswith(" / build"):
            continue

        key = None
        if "windows-latest-msvc" in name:
            key = "msvc"
        elif "ubuntu-latest" in name:
            key = "ubuntu"
        elif "windows-latest-rtools" in name:
            key = "rtools"
        elif "macos-latest" in name:
            key = "macos"

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
    """Generates 4 stacked subplots (one per platform) with independent Y-scales."""
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

    platforms = [
        ("MSVC", [to_min(d["msvc"]) for d in history_data], "#1f77b4"),
        ("Ubuntu", [to_min(d["ubuntu"]) for d in history_data], "#ff7f0e"),
        ("Rtools", [to_min(d["rtools"]) for d in history_data], "#2ca02c"),
        ("macOS", [to_min(d["macos"]) for d in history_data], "#d62728"),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    title_suffix = f" (Test: {test_name})" if test_name else " (Job Total)"
    fig.suptitle(
        f"gstlearn - Execution Times{title_suffix}", fontsize=14, fontweight="bold"
    )

    for i, (name, vals, color) in enumerate(platforms):
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

    filename = f"nonreg_times_{test_name}.png" if test_name else "nonreg_times.png"
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

    headers = ["Run", "Date", "MSVC", "Ubuntu", "Rtools", "macOS"]
    widths = [8, 13, 14, 14, 14, 14]

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
            msvc_sec = test_times.get("msvc")
            ubuntu_sec = test_times.get("ubuntu")
            rtools_sec = test_times.get("rtools")
            macos_sec = test_times.get("macos")

            msvc_val = f"{msvc_sec:.2f}s" if msvc_sec is not None else "-"
            ubuntu_val = f"{ubuntu_sec:.2f}s" if ubuntu_sec is not None else "-"
            rtools_val = f"{rtools_sec:.2f}s" if rtools_sec is not None else "-"
            macos_val = f"{macos_sec:.2f}s" if macos_sec is not None else "-"
        else:
            jobs = get_build_jobs(internal_id)
            msvc_sec = jobs.get("msvc", {}).get("elapsed")
            ubuntu_sec = jobs.get("ubuntu", {}).get("elapsed")
            rtools_sec = jobs.get("rtools", {}).get("elapsed")
            macos_sec = jobs.get("macos", {}).get("elapsed")

            msvc_val = format_elapsed(msvc_sec)
            ubuntu_val = format_elapsed(ubuntu_sec)
            rtools_val = format_elapsed(rtools_sec)
            macos_val = format_elapsed(macos_sec)

        history_data.append(
            {
                "run": run["run_number"],
                "date": run_date,
                "msvc": msvc_sec,
                "ubuntu": ubuntu_sec,
                "rtools": rtools_sec,
                "macos": macos_sec,
            }
        )

        values = [
            str(run["run_number"]),
            run_date,
            msvc_val,
            ubuntu_val,
            rtools_val,
            macos_val,
        ]

        print("".join(f"{v:<{w}}" for v, w in zip(values, widths)), flush=True)

    print()
    generate_chart(history_data, test_name)


if __name__ == "__main__":
    main()
