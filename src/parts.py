#!/usr/bin/env python3
"""parts.py — the part list for the EnzymeSwitch device (Module 1).

A "part" is one piece of DNA (a promoter, a gene, ...).
A "device" is those parts snapped together in a specific order.
This file stores the list of parts and prints it when you run it.
"""

from dataclasses import dataclass

# ----------------------------------------------------------------------
# THE DESIGN — the one true order of our device.
# Tuple (not list) = cannot be accidentally changed later. Safety choice.
# ----------------------------------------------------------------------
DEVICE_ORDER = (
    "P_con_lacI",   # 1 - always-ON promoter that powers the repressor
    "P_t7_lacO",    # 2 - THE main switch (OFF until IPTG arrives)
    "RBS",          # 3 - ribosome grab-handle (starts protein making)
    "CDS_PETase",   # 4 - the enzyme gene, our output
    "T_terminator", # 5 - end-of-line signal for transcription
    "Backbone",     # 6 - the plasmid chassis that holds everything
)

# ----------------------------------------------------------------------
# A Part is one DNA piece. @dataclass writes the boring boilerplate for
# us (constructor, comparison, printing) — we just list the fields.
# ----------------------------------------------------------------------
@dataclass
class Part:
    id: str                     # e.g. "BBa_R0010" — ID in a registry
    part_type: str = "unknown"  # promoter / RBS / CDS / terminator / backbone
    name: str = ""              # human-readable name
    bb_id: str = ""             # BioBrick / registry ID, e.g. "BBa_J23100"    
    seq: str | None = None      # DNA letters. LEFT EMPTY for now. Gonna get from databanks.
    source: str = ""            # where the part comes from (URL / GenBank)
    role: str = ""              # what this part does in OUR device

# ----------------------------------------------------------------------
# THE ALARM — checks that a list of parts follows DEVICE_ORDER.
# ----------------------------------------------------------------------
def validate_device_order(parts):
    """Check that part ids follow DEVICE_ORDER. Raise ValueError if not."""
    # pull just the ids into a plain list for easy comparison
    actual = [p.id for p in parts]

    # the heart of the check
    if actual != list(DEVICE_ORDER):
        # find which slot differs FIRST, so the message is useful
        for i, (got, want) in enumerate(zip(actual, DEVICE_ORDER), start=1):
            if got != want:
                break
        raise ValueError(
            f"Device order mismatch at slot {i}: got '{got}' but expected '{want}'.\n"
            f"   actual list:  {actual}\n"
            f"   design wants: {list(DEVICE_ORDER)}"
        )

    # survived to here = order is perfect
    print(
        f"Device order assertion: PASS "
        f"({len(parts)}/{len(DEVICE_ORDER)} slots in canonical order)"
    )

# ----------------------------------------------------------------------
# THE BUILDER — makes the 6 real parts and runs the alarm on them,
# so every future module gets a CHECKED device, never a broken one.
# ----------------------------------------------------------------------
def device_parts():
    """Build the 6 parts of the EnzymeSwitch device and return them in order."""
    parts = [
        Part(
            id="P_con_lacI",    # logical slot name (this is what DEVICE_ORDER checks)
            bb_id="BBa_J23100", # real iGEM part to order
            part_type="promoter",
            name="constitutive promoter (lacI)",
            source="iGEM registry — VERIFY at build time",
            role="always-ON: keeps the lacI repressor produced",
        ),
        Part(
            id="P_t7_lacO",
            bb_id="BBa_R0010",
            part_type="promoter",
            name="T7 promoter + lacO",
            source="iGEM registry — VERIFY at build time",
            role="main switch: OFF until IPTG arrives",
        ),
        Part(
            id="RBS",
            bb_id="BBa_B0034",
            part_type="RBS",
            name="strong RBS",
            source="iGEM registry — VERIFY at build time",
            role="translation initiation: ribosome grab-handle",
        ),
        Part(
            id="CDS_PETase",
            bb_id="PETase_Ec",
            part_type="CDS",
            name="PETase CDS (Ideonella sakaiensis)",
            source="UniProt A0A0K8P6T7 — VERIFY at fetch (Module 4)",
            role="output: the plastic-degrading enzyme",
        ),
        Part(
            id="T_terminator",
            bb_id="BBa_B0015",
            part_type="terminator",
            name="double terminator",
            source="iGEM registry — VERIFY at build time",
            role="stops transcription cleanly",
        ),
        Part(
            id="Backbone",
            bb_id="pET-28a(+)",
            part_type="backbone",
            name="pET-28a(+) expression vector",
            source="Novagen — sequence from GenBank at build time",
            role="chassis: replication origin + KanR + lacI gene",
        ),
    ]
    # THE GATE — run the alarm; if order is broken, this raises and stops.
    validate_device_order(parts)

    return parts

# ----------------------------------------------------------------------
# Runs only when you type:  python src/parts.py
# (not when another file imports this one)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parts = device_parts()

    print("\n=== EnzymeSwitch device map (Module 1) ===\n")
    # headers, padded so the table lines up
    print(f"{'#':<3}{'Slot':<14}{'Part ID':<14}{'Name':<30}{'Type':<12}Role")
    print("-" * 88)

    for i, part in enumerate(parts, start=1):
        print(
            f"{i:<3}{part.id:<14}{part.name:<30}"
            f"{part.part_type:<12}{part.role}"
        )

    # honest line: no real sequences yet, and we say so out loud
    print("\nSequences: NOT fetched yet (integrity policy — fetch at build time).")
