#!/usr/bin/env python
# Author: The oslumen Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at doc.oslu.men/pyddscat

import importlib.metadata as metadata


def get_meta(metadata):
    try:
        data = metadata.metadata("pyddscat")
        __version__ = metadata.version("pyddscat")
        __author__ = data.get("author")
        __description__ = data.get("summary")
    except Exception:
        data = dict(License="unknown")
        __version__ = "unknown"
        __author__ = "unknown"
        __description__ = "unknown"
    return __version__, __author__, __description__, data


__version__, __author__, __description__, data = get_meta(metadata)
