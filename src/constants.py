"""Local PDG decay table.

pdgData is a dict keyed by PDG id. Each record has "name", "charge",
"mass" (MeV), "decay", and "lifetime". A decay entry is
[name, daughter PDG ids..., branching fraction, cumulative].
Example: pdgData[211]["decay"] begins
[["munu", -13, 14, 0.999877, 0.999877], ...].

Unknown PDG ids raise UnknownPDGId.

The records are the 2016-07-09 Wayback Machine snapshot of
http://frankliuao.com/downloads/Codes/pdg_py.dict, stored in
pdg_table.json next to this file. Nothing is downloaded and nothing
is passed to exec.
"""

import json
import os

try:
    _STRING_TYPES = (str, unicode)
except NameError:
    _STRING_TYPES = (str,)

_TABLE_FILENAME = "pdg_table.json"


class UnknownPDGId(KeyError):
    """Raised when a PDG id is not in the vendored decay table."""


def _coerce_pdg_id(pdg_id):
    """Return an int PDG id, accepting integral floats such as 211.0.

    Beam arrays store PDG ids as floats. Non-integral values and
    non-numeric values are unknown ids.
    """
    if isinstance(pdg_id, bool) or isinstance(pdg_id, _STRING_TYPES):
        raise UnknownPDGId("Unknown PDG id: %r" % (pdg_id,))
    try:
        number = float(pdg_id)
        key = int(number)
    except (TypeError, ValueError, OverflowError):
        raise UnknownPDGId("Unknown PDG id: %r" % (pdg_id,))
    if number != key:
        raise UnknownPDGId("Unknown PDG id: %r" % (pdg_id,))
    return key


class PDGTable(dict):
    """PDG id -> particle record. Missing ids raise UnknownPDGId."""

    def __getitem__(self, pdg_id):
        key = _coerce_pdg_id(pdg_id)
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            known = ", ".join(str(item) for item in sorted(self.keys()))
            raise UnknownPDGId(
                "Unknown PDG id: %s. It is not in the vendored decay table. "
                "Known ids: %s" % (key, known)
            )


def _default_table_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        _TABLE_FILENAME)


def load_pdg_data(path=None):
    """Load the vendored PDG table from disk.

    Returns (pdg_data, pdg_map). pdg_data raises UnknownPDGId for ids
    that are not in the file. The file is not fetched from the network.
    """
    path = _default_table_path() if path is None else path
    try:
        handle = open(path, "r")
    except IOError:
        raise IOError(
            "Vendored PDG table not found: %s. "
            "The decay table is local data and is not downloaded." % path
        )
    try:
        try:
            payload = json.load(handle)
        except ValueError as exc:
            raise ValueError(
                "Could not parse vendored PDG table %s: %s" % (path, exc)
            )
    finally:
        handle.close()

    if not isinstance(payload, dict) or "pdgData" not in payload \
            or "pdgMap" not in payload:
        raise ValueError(
            "Vendored PDG table %s is missing pdgData or pdgMap." % path
        )
    raw_data = payload["pdgData"]
    raw_map = payload["pdgMap"]
    if not isinstance(raw_data, dict) or len(raw_data) == 0:
        raise ValueError(
            "Vendored PDG table %s has an empty pdgData object." % path
        )
    if not isinstance(raw_map, dict):
        raise ValueError(
            "Vendored PDG table %s has a pdgMap that is not an object." % path
        )

    table = PDGTable()
    for key, record in raw_data.items():
        if not isinstance(key, _STRING_TYPES):
            raise ValueError(
                "Vendored PDG table has a non-integer id: %r" % (key,)
            )
        try:
            pdg_id = int(key)
        except (TypeError, ValueError):
            raise ValueError(
                "Vendored PDG table has a non-integer id: %r" % (key,)
            )
        if str(pdg_id) != str(key):
            raise ValueError(
                "Vendored PDG table has a non-integer id: %r" % (key,)
            )
        if not isinstance(record, dict):
            raise ValueError(
                "Vendored PDG record for %s is not an object." % pdg_id
            )
        for field in ("name", "charge", "mass", "decay"):
            if field not in record:
                raise ValueError(
                    "Vendored PDG record for %s is missing %s."
                    % (pdg_id, field)
                )
        table[pdg_id] = record
    return table, raw_map


pdgData, pdgMap = load_pdg_data()
