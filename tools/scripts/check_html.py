import requests
from bs4 import BeautifulSoup
from urllib.parse import unquote
import os

COURSES_PYTHON_FOLDER = (
    "https://soft.mines-paristech.fr/gstlearn/courses-latest/python/"
)
COURSES_R_FOLDER = "https://soft.mines-paristech.fr/gstlearn/courses-latest/r/"
COURSES_PAGE = "https://gstlearn.org/?page_id=193"

DEMOS_PYTHON_FOLDER = "https://soft.mines-paristech.fr/gstlearn/demos-latest/python/"
DEMOS_R_FOLDER = "https://soft.mines-paristech.fr/gstlearn/demos-latest/r/"
DEMOS_PYTHON_PAGE = "https://gstlearn.org/?page_id=187"
DEMOS_R_PAGE = "https://gstlearn.org/?page_id=420"

DEMOS_BLACK_LIST = {
    "Tuto_Block_Variance.html",
    "Tuto_Animate.html",
    "Tuto_Kriging.html",
    "Tuto_Model.html",
    "Tuto_Variogram.html",
    "Tuto_Widget.html",
    "Inheritance.html",
    "Starting.html",
}

session = requests.Session()


def get_soup(url):
    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
        return BeautifulSoup(response.text, "html.parser")
    except requests.exceptions.RequestException as e:
        print(f"Connection error on {url} : {e}")
        return None


def extract_html_filenames(url, condition_func):
    soup = get_soup(url)
    if not soup:
        return set()

    files = set()
    for link in soup.find_all("a"):
        href = link.get("href")
        if href and condition_func(href):
            filename = unquote(href.split("/")[-1])
            files.add(filename)
    return files


def extract_html_links_of_file(url, condition_func):
    soup = get_soup(url)
    if not soup:
        return set()

    links = set()
    for link in soup.find_all("a"):
        href = link.get("href")
        if href and condition_func(href):
            decoded_href = unquote(href)
            links.add(decoded_href)

    return links


def check_courses(python_folder_url, r_folder_url, shared_page_url):
    print(f"\n{'=' * 50}")
    print("-- Analysis report (Courses)")
    print(f"{'=' * 50}")

    python_files = extract_html_filenames(
        python_folder_url, lambda h: h.endswith(".html")
    )
    r_files = extract_html_filenames(r_folder_url, lambda h: h.endswith(".html"))
    all_page_links = extract_html_filenames(shared_page_url, lambda h: ".html" in h)

    linked_files = extract_html_links_of_file(
        shared_page_url, lambda h: h.endswith(".html")
    )

    wrong_folder = False
    for link in linked_files:
        copy_link = link

        folder, file = os.path.split(copy_link)

        folder_clean = folder.strip("/")

        if folder_clean == COURSES_PYTHON_FOLDER.strip(
            "/"
        ) or folder_clean == COURSES_R_FOLDER.strip("/"):
            continue

        print(f"Wrong folder for file {file}")
        wrong_folder = True

    if not python_files and not r_files and not all_page_links:
        return

    missing_python = python_files - all_page_links
    missing_r = r_files - all_page_links
    all_folder_files = python_files | r_files
    dead_links = all_page_links - all_folder_files

    if not missing_python and not missing_r and not dead_links and not wrong_folder:
        print("Everything is perfectly in sync!")
        return

    if missing_python:
        print("\nPython (In the folder but not on the page) :")
        for f in sorted(missing_python):
            print(f"  [+] {f}")

    if missing_r:
        print("\nR (In the folder but not on the page) :")
        for f in sorted(missing_r):
            print(f"  [+] {f}")

    if dead_links:
        print("\nTo be corrected (Broken links on the page) :")
        for f in sorted(dead_links):
            print(f"  [-] {f}")


def check_demos(folder_url, page_url, language):
    print(f"\n{'=' * 50}")
    print(f"-- Analysis report (Demos {language})")
    print(f"{'=' * 50}")

    folder_files = extract_html_filenames(folder_url, lambda h: h.endswith(".html"))
    page_links = extract_html_filenames(page_url, lambda h: ".html" in h)

    linked_files = extract_html_links_of_file(page_url, lambda h: h.endswith(".html"))

    wrong_folder = False
    for link in linked_files:
        copy_link = link

        folder, file = os.path.split(copy_link)

        folder_clean = folder.strip("/")

        if folder_clean == DEMOS_PYTHON_FOLDER.strip(
            "/"
        ) or folder_clean == DEMOS_R_FOLDER.strip("/"):
            continue

        print(f"Wrong folder for file {file}")
        wrong_folder = True

    if not folder_files and not page_links:
        return

    missing_on_page = folder_files - page_links - DEMOS_BLACK_LIST
    dead_links = page_links - folder_files

    if not missing_on_page and not dead_links and not wrong_folder:
        print("Everything is perfectly in sync!")
        return

    if missing_on_page:
        print("\nIn the folder but not on the page :")
        for f in sorted(missing_on_page):
            print(f"  [+] {f}")

    if dead_links:
        print("\nTo be corrected (Broken links on the page) :")
        for f in sorted(dead_links):
            print(f"  [-] {f}")


if __name__ == "__main__":
    check_courses(COURSES_PYTHON_FOLDER, COURSES_R_FOLDER, COURSES_PAGE)

    check_demos(DEMOS_PYTHON_FOLDER, DEMOS_PYTHON_PAGE, "Python")
    check_demos(DEMOS_R_FOLDER, DEMOS_R_PAGE, "R")

    print(f"\n{'=' * 50}")
    print("-- Analysis complete.")
    print(f"{'=' * 50}")
