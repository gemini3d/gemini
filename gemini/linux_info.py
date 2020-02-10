#!/usr/bin/env python
"""
get Linux system info
"""
from configparser import ConfigParser
import typing as T
from pathlib import Path


def os_release() -> T.List[str]:
    """
    reads /etc/os-release with fallback to legacy methods

    returns
    -------
    'rhel' or 'debian'
    """
    fn = Path("/etc/os-release")
    if not fn.is_file():
        if Path("/etc/redhat-release").is_file() or Path("/etc/centos-release").is_file():
            return ["rhel"]
        elif Path("/etc/debian_version").is_file():
            return ["debian"]

    C = ConfigParser(inline_comment_prefixes=("#", ";"))
    ini = "[all]" + fn.read_text()
    C.read_string(ini)
    return C["all"].get("ID_LIKE").strip('"').strip("'").split()


def get_package_manager(like: T.List[str] = None) -> str:
    if not like:
        like = os_release()
    if isinstance(like, str):
        like = [like]

    if {"centos", "rhel", "fedora"}.intersection(like):
        return "yum"
    elif {"debian", "ubuntu"}.intersection(like):
        return "apt"
    else:
        raise ValueError(f"Unknown ID_LIKE={like}, please file bug report or manually specify package manager")
