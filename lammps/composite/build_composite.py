#!/usr/bin/env python3
"""Build a Fe3O4/UiO-66/Fe3O4 composite (Section 3.2 of the article).

Steps
-----
1. Read a relaxed UiO-66 cell (LAMMPS data file, atom_style full, with the
   UFF topology and DDEC charges), replicated in z as needed.
2. Cut a slab normal to z. The cut is made by building units, not by atoms:
   a Zr6O4(OH)4 node is kept if all of its atoms lie inside [zlo, zhi], and a
   BDC linker is kept only if all of its atoms lie inside the window and both
   nodes to which it is coordinated are kept. In this way only Zr-O(carboxylate)
   bonds are broken by the cut.
3. Saturate every Zr that lost a carboxylate oxygen with a hydroxyl group
   (Zr-OH cap). The cap oxygen is placed along the original Zr-O bond, and the
   topology of the caps (bonds, angles, and dihedrals) is assigned by analogy
   with the bulk framework.
4. Place two Fe3O4 slabs symmetrically below and above the UiO-66 slab, with a
   gap between the surfaces, and add vacuum at both ends of the box.
5. Write the composite as a LAMMPS data file (atom_style full). The atom types
   of Fe3O4 are appended after those of the framework.

Usage
-----
    python build_composite.py --mof {mof.data} --oxide {Fe3O4_slab.data} \
        --zlo {zlo} --zhi {zhi} --gap {gap} --vacuum {vacuum} --out {composite.data}

All lengths are in Angstrom. Run with -h for all options. The script prints
an audit (atom counts, net charge, shortest contacts) that should be checked
before the composite is used.
"""
import argparse
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

SECTIONS_COEFF = ["Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
                  "Dihedral Coeffs", "Improper Coeffs"]
SECTIONS_TOPO = {"Bonds": 2, "Angles": 3, "Dihedrals": 4, "Impropers": 4}


# --------------------------------------------------------------------------
# LAMMPS data file I/O (orthogonal box, atom_style full)
# --------------------------------------------------------------------------
def read_data(path):
    lines = open(path).read().splitlines()
    d = {"box": {}, "coeff": {}, "counts": {}, "topo": {k: [] for k in SECTIONS_TOPO}}
    sec, i = None, 0
    header_done = False
    for raw in lines[1:]:
        s = raw.split("#")[0].strip()
        if not s:
            continue
        m = re.match(r"^([-\d.eE+]+)\s+([-\d.eE+]+)\s+([xyz])lo\s+\3hi$", s)
        if m:
            d["box"][m.group(3)] = (float(m.group(1)), float(m.group(2)))
            continue
        m = re.match(r"^(\d+)\s+(atoms|bonds|angles|dihedrals|impropers|atom types|bond types|angle types|dihedral types|improper types)$", s)
        if m and not header_done:
            d["counts"][m.group(2)] = int(m.group(1))
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", s):
            sec = s
            header_done = True
            continue
        if sec in SECTIONS_COEFF:
            d["coeff"].setdefault(sec, []).append(raw.strip())
        elif sec == "Atoms":
            p = s.split()
            d.setdefault("atoms", []).append(
                (int(p[0]), int(p[1]), int(p[2]), float(p[3]),
                 float(p[4]), float(p[5]), float(p[6])))
        elif sec in SECTIONS_TOPO:
            p = list(map(int, s.split()))
            d["topo"][sec].append((p[1], tuple(p[2:2 + SECTIONS_TOPO[sec]])))
    a = sorted(d["atoms"])
    d["id"] = np.array([r[0] for r in a])
    d["mol"] = np.array([r[1] for r in a])
    d["type"] = np.array([r[2] for r in a])
    d["q"] = np.array([r[3] for r in a])
    d["xyz"] = np.array([r[4:7] for r in a])
    d["L"] = np.array([d["box"][k][1] - d["box"][k][0] for k in "xyz"])
    d["lo"] = np.array([d["box"][k][0] for k in "xyz"])
    d["xyz"] = d["xyz"] - d["lo"]
    d["xyz"] = d["xyz"] - np.floor(d["xyz"] / d["L"]) * d["L"]
    return d


def write_data(path, title, L, types, mol, q, xyz, coeff, topo, ntypes):
    with open(path, "w") as f:
        f.write(title + "\n\n")
        f.write(f"{len(types)} atoms\n")
        for k, name in [("Bonds", "bonds"), ("Angles", "angles"),
                        ("Dihedrals", "dihedrals"), ("Impropers", "impropers")]:
            f.write(f"{len(topo[k])} {name}\n")
        f.write("\n")
        for k in ["atom", "bond", "angle", "dihedral", "improper"]:
            f.write(f"{ntypes[k]} {k} types\n")
        f.write("\n")
        for k, n in zip("xyz", L):
            f.write(f"0.0 {n:.6f} {k}lo {k}hi\n")
        for sec in SECTIONS_COEFF:
            if coeff.get(sec):
                hdr = sec + (" # lj/cut/coul/long" if sec == "Pair Coeffs" else "")
                f.write(f"\n{hdr}\n\n" + "\n".join(coeff[sec]) + "\n")
        f.write("\nAtoms # full\n\n")
        for i in range(len(types)):
            f.write(f"{i + 1} {mol[i]} {types[i]} {q[i]:.6f} "
                    f"{xyz[i, 0]:.6f} {xyz[i, 1]:.6f} {xyz[i, 2]:.6f}\n")
        for sec in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
            if topo[sec]:
                f.write(f"\n{sec}\n\n")
                for j, (t, ids) in enumerate(topo[sec], 1):
                    f.write(f"{j} {t} " + " ".join(str(x + 1) for x in ids) + "\n")


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def mic(v, L, periodic):
    v = v.copy()
    for k in range(3):
        if periodic[k]:
            v[..., k] -= L[k] * np.round(v[..., k] / L[k])
    return v


def angle_deg(a, b, c, xyz, L, per):
    u = mic(xyz[a] - xyz[b], L, per)
    w = mic(xyz[c] - xyz[b], L, per)
    cs = np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w)
    return np.degrees(np.arccos(np.clip(cs, -1, 1)))


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mof", required=True, help="relaxed UiO-66 cell (atom_style full, UFF topology, DDEC charges)")
    ap.add_argument("--oxide", required=True, help="relaxed Fe3O4 slab (atom_style full)")
    ap.add_argument("--zlo", type=float, required=True, help="lower bound of the slab window (A, from the bottom of the cell)")
    ap.add_argument("--zhi", type=float, required=True, help="upper bound of the slab window (A)")
    ap.add_argument("--gap", type=float, required=True, help="gap between the MOF and oxide surfaces (A)")
    ap.add_argument("--vacuum", type=float, required=True, help="vacuum at each end of the box (A)")
    ap.add_argument("--out", required=True, help="output data file")
    ap.add_argument("--types", default="1,2,3,4,5",
                    help="framework types C_R,H_,O_R,O_3_z,Zr (default 1,2,3,4,5, as in data.UiO-66.ddec)")
    ap.add_argument("--r_zro", type=float, default=2.15, help="Zr-O(cap) distance (A)")
    ap.add_argument("--r_oh", type=float, default=0.97, help="O-H(cap) distance (A)")
    ap.add_argument("--ang_zroh", type=float, default=110.0, help="Zr-O-H(cap) angle (degrees)")
    ap.add_argument("--neutralize", choices=["caps", "none"], default="caps",
                    help="spread the residual charge over the cap atoms (default) or leave it")
    ap.add_argument("--lateral_tol", type=float, default=0.02,
                    help="largest accepted relative mismatch of the lateral box lengths of MOF and oxide")
    a = ap.parse_args()

    tC, tH, tOR, tO3, tZr = map(int, a.types.split(","))
    m = read_data(a.mof)
    L, xyz, typ, q = m["L"], m["xyz"], m["type"], m["q"]
    N = len(typ)
    idx = {aid: i for i, aid in enumerate(m["id"])}
    per = (True, True, True)

    # ---- bond graph and building units -------------------------------------
    nbr = defaultdict(set)
    bond_type = {}
    for t, (i, j) in m["topo"]["Bonds"]:
        i, j = idx[i], idx[j]
        nbr[i].add(j); nbr[j].add(i)
        bond_type[(typ[i], typ[j])] = t
        bond_type[(typ[j], typ[i])] = t
    zr_or = lambda i, j: {typ[i], typ[j]} == {tZr, tOR}
    unit = -np.ones(N, int)
    nunit = 0
    for s in range(N):
        if unit[s] >= 0:
            continue
        stack = [s]; unit[s] = nunit
        while stack:
            i = stack.pop()
            for j in nbr[i]:
                if unit[j] < 0 and not zr_or(i, j):
                    unit[j] = nunit; stack.append(j)
        nunit += 1
    members = defaultdict(list)
    for i, u in enumerate(unit):
        members[u].append(i)
    is_node = {u: any(typ[i] == tZr for i in at) for u, at in members.items()}
    print(f"building units: {sum(is_node.values())} nodes, {nunit - sum(is_node.values())} linkers")

    # unwrap each unit around its first atom
    unwrapped = xyz.copy()
    for u, at in members.items():
        root = at[0]; seen = {root}; stack = [root]
        while stack:
            i = stack.pop()
            for j in nbr[i]:
                if unit[j] == u and j not in seen:
                    unwrapped[j] = unwrapped[i] + mic(xyz[j] - xyz[i], L, per)
                    seen.add(j); stack.append(j)
        c = unwrapped[at].mean(axis=0)
        shift = np.floor(c / L) * L
        unwrapped[at] -= shift

    inside = {u: bool(np.all((unwrapped[at, 2] >= a.zlo) & (unwrapped[at, 2] <= a.zhi)))
              for u, at in members.items()}
    keep_unit = {u for u in members if is_node[u] and inside[u]}
    for u, at in members.items():
        if is_node[u] or not inside[u]:
            continue
        partners = {unit[j] for i in at for j in nbr[i] if zr_or(i, j)}
        if partners and partners <= keep_unit:
            keep_unit.add(u)
    keep = np.array(sorted(i for u in keep_unit for i in members[u]))
    if len(keep) == 0:
        sys.exit("empty slab: widen [zlo, zhi]")
    print(f"kept {len(keep)} of {N} atoms "
          f"({sum(is_node[u] for u in keep_unit)} nodes, {sum(not is_node[u] for u in keep_unit)} linkers)")

    # ---- caps ---------------------------------------------------------------
    kept = set(keep.tolist())
    pos = unwrapped
    caps = []  # (zr_index, O_xyz, H_xyz)
    for i in keep:
        if typ[i] != tZr:
            continue
        cen = pos[members[unit[i]]].mean(axis=0)
        for j in nbr[i]:
            if typ[j] == tOR and j not in kept:
                u = mic(xyz[j] - xyz[i], L, per); u /= np.linalg.norm(u)
                o = pos[i] + a.r_zro * u
                w = (o - cen) - np.dot(o - cen, u) * u
                if np.linalg.norm(w) < 1e-6:
                    w = np.cross(u, [0, 0, 1.0]) if abs(u[2]) < 0.9 else np.cross(u, [1.0, 0, 0])
                w /= np.linalg.norm(w)
                th = np.radians(180.0 - a.ang_zroh)
                h = o + a.r_oh * (np.cos(th) * u + np.sin(th) * w)
                caps.append((i, o, h))
    print(f"Zr-OH caps added: {len(caps)}")

    # new atom arrays (MOF part)
    old2new = {old: k for k, old in enumerate(keep)}
    m_types = list(typ[keep]); m_q = list(q[keep]); m_xyz = [pos[i] for i in keep]; m_mol = list(m["mol"][keep])
    qO = q[typ == tO3].mean()
    qH = np.mean([q[i] for i in range(N) if typ[i] == tH and any(typ[j] == tO3 for j in nbr[i])] or [q[typ == tH].mean()])
    cap_atoms = []
    new_bonds = []
    for zr, o, h in caps:
        io = len(m_types); m_types.append(tO3); m_q.append(qO); m_xyz.append(o); m_mol.append(m["mol"][zr])
        ih = len(m_types); m_types.append(tH); m_q.append(qH); m_xyz.append(h); m_mol.append(m["mol"][zr])
        cap_atoms += [io, ih]
        new_bonds += [(old2new[zr], io), (io, ih)]
    m_types = np.array(m_types); m_q = np.array(m_q); m_xyz = np.array(m_xyz); m_mol = np.array(m_mol)

    # ---- topology -----------------------------------------------------------
    topo = {k: [] for k in SECTIONS_TOPO}
    for sec in SECTIONS_TOPO:
        for t, ids in m["topo"][sec]:
            ii = [idx[x] for x in ids]
            if all(x in kept for x in ii):
                topo[sec].append((t, tuple(old2new[x] for x in ii)))
    per_slab = (True, True, False)
    # bonds of the caps
    for i, j in new_bonds:
        t = bond_type.get((m_types[i], m_types[j]))
        if t is None:
            sys.exit(f"no bond type for types {m_types[i]}-{m_types[j]}")
        topo["Bonds"].append((t, (i, j)))
    # angle types by analogy: same end and centre types, closest mean angle
    ang_ref = defaultdict(list)
    for t, (i, j, k) in m["topo"]["Angles"]:
        i, j, k = idx[i], idx[j], idx[k]
        key = (min(typ[i], typ[k]), typ[j], max(typ[i], typ[k]))
        ang_ref[key].append((t, angle_deg(i, j, k, xyz, L, per)))
    ang_mean = {key: {t: np.mean([x for tt, x in v if tt == t]) for t in {tt for tt, _ in v}} for key, v in ang_ref.items()}
    adj = defaultdict(set)
    for _, (i, j) in topo["Bonds"]:
        adj[i].add(j); adj[j].add(i)
    capset = set(cap_atoms)
    nmiss = 0
    for c in set(x for b in new_bonds for x in b):
        nb = sorted(adj[c])
        for p in range(len(nb)):
            for r in range(p + 1, len(nb)):
                i, k = nb[p], nb[r]
                if not ({i, c, k} & capset):
                    continue
                key = (min(m_types[i], m_types[k]), m_types[c], max(m_types[i], m_types[k]))
                if key not in ang_mean:
                    nmiss += 1; continue
                th = angle_deg(i, c, k, m_xyz, L, per_slab)
                t = min(ang_mean[key], key=lambda tt: abs(ang_mean[key][tt] - th))
                topo["Angles"].append((t, (i, c, k)))
    # dihedral types by analogy (type quadruplet)
    dih_ref = {}
    for t, ids in m["topo"]["Dihedrals"]:
        k4 = tuple(typ[idx[x]] for x in ids)
        dih_ref.setdefault(k4, t); dih_ref.setdefault(k4[::-1], t)
    seen = set()
    for b, c in [e for e in new_bonds] + [(j, i) for i, j in new_bonds]:
        for aa in adj[b]:
            for dd in adj[c]:
                if len({aa, b, c, dd}) < 4:
                    continue
                key = (aa, b, c, dd) if aa < dd else (dd, c, b, aa)
                if key in seen:
                    continue
                seen.add(key)
                t = dih_ref.get(tuple(m_types[x] for x in key))
                if t is not None:
                    topo["Dihedrals"].append((t, key))
    if nmiss:
        print(f"note: {nmiss} cap angles have no analogue in the bulk and were omitted")

    # ---- charge -------------------------------------------------------------
    net = m_q.sum()
    print(f"net charge of the capped slab: {net:+.4f} e")
    if a.neutralize == "caps" and cap_atoms:
        m_q[cap_atoms] -= net / len(cap_atoms)
        print(f"residual spread over {len(cap_atoms)} cap atoms; net charge now {m_q.sum():+.2e} e")

    # ---- oxide slabs --------------------------------------------------------
    ox = read_data(a.oxide)
    mis = (ox["L"][:2] - L[:2]) / L[:2]
    print(f"lateral mismatch oxide/MOF: x {100 * mis[0]:+.2f} %, y {100 * mis[1]:+.2f} %")
    if np.any(np.abs(mis) > a.lateral_tol):
        sys.exit("lateral mismatch above --lateral_tol: replicate the Fe3O4 cell so that its lateral "
                 "box lengths match those of the MOF slab")
    oxyz = ox["xyz"].copy()
    oxyz[:, :2] = oxyz[:, :2] / ox["L"][:2] * L[:2]
    oxyz[:, 2] -= oxyz[:, 2].min()
    t_ox = oxyz[:, 2].max()
    nmof_t = m["counts"]["atom types"]

    m_xyz[:, 2] -= m_xyz[:, 2].min()
    t_mof = m_xyz[:, 2].max()
    z_bot = a.vacuum
    z_mof = z_bot + t_ox + a.gap
    z_top = z_mof + t_mof + a.gap
    Lz = z_top + t_ox + a.vacuum
    m_xyz[:, 2] += z_mof
    m_xyz[:, :2] -= np.floor(m_xyz[:, :2] / L[:2]) * L[:2]
    bot = oxyz + [0, 0, z_bot]
    top = oxyz + [0, 0, z_top]

    nox = len(ox["type"])
    maxmol = int(m_mol.max())
    types = np.concatenate([m_types, ox["type"] + nmof_t, ox["type"] + nmof_t])
    charges = np.concatenate([m_q, ox["q"], ox["q"]])
    coords = np.vstack([m_xyz, bot, top])
    mols = np.concatenate([m_mol, np.full(nox, maxmol + 1), np.full(nox, maxmol + 2)])

    coeff = {k: list(v) for k, v in m["coeff"].items()}
    for sec in ["Masses", "Pair Coeffs"]:
        for line in ox["coeff"].get(sec, []):
            p = line.split("#")[0].split()
            coeff.setdefault(sec, []).append(" ".join([str(int(p[0]) + nmof_t)] + p[1:]))
    ntypes = {"atom": nmof_t + ox["counts"]["atom types"],
              "bond": m["counts"]["bond types"], "angle": m["counts"]["angle types"],
              "dihedral": m["counts"].get("dihedral types", 0), "improper": m["counts"].get("improper types", 0)}
    Lout = np.array([L[0], L[1], Lz])

    # ---- audit --------------------------------------------------------------
    tree = cKDTree(coords, boxsize=[L[0], L[1], Lz * 10])
    pairs = tree.query_pairs(1.0)
    bonded = {(min(i, j), max(i, j)) for _, (i, j) in topo["Bonds"]}
    close = [p for p in pairs if p not in bonded]
    nm = len(m_types)
    d_if, _ = cKDTree(coords[nm:], boxsize=[L[0], L[1], Lz * 10]).query(coords[:nm])
    print(f"composite: {len(types)} atoms ({nm} MOF, {2 * nox} Fe3O4), box {L[0]:.3f} x {L[1]:.3f} x {Lz:.3f} A")
    print(f"MOF slab thickness {t_mof:.2f} A; oxide slab thickness {t_ox:.2f} A")
    print(f"shortest MOF-oxide distance {d_if.min():.2f} A; nonbonded pairs closer than 1.0 A: {len(close)}")
    print(f"net charge of the composite: {charges.sum():+.4f} e")
    mass = {int(l.split()[0]): float(l.split()[1]) for l in coeff["Masses"]}
    m_mof = sum(mass[t] for t in m_types); m_ox = sum(mass[t] for t in types[nm:])
    print(f"Fe3O4 weight fraction: {100 * m_ox / (m_ox + m_mof):.1f} wt %")

    write_data(a.out, "Fe3O4/UiO-66/Fe3O4 composite built with build_composite.py",
               Lout, types, mols, charges, coords, coeff, topo, ntypes)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
