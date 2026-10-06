"""Independent evaluator of the MultiSurf deck Input/C1.ms2 (test-only, no MATLAB code involved).

Reads the deck text and evaluates, for given heights z, the horizontal section of the hull:
area, second moments Ixx = integral y^2 dA and Iyy = integral x^2 dA about the body origin
(x = y = 0), and the x half-width.

Deck semantics used (read from the deck text):
  FramePoint   point = referenced point + (dx, dy, dz); '*' is the origin.
  MirrPoint    mirror image of a point in the stated plane (X=0).
  BCurve       clamped uniform B-spline of the stated degree through its control points,
               parameter t in [0, 1].
  Arc          '/ * 2 a c b': circular arc from a to b about centre c (the arc of angle <= 180 deg).
  PolyCurve2   concatenation of the listed curves, each taking an equal share of t in [0, 1].
  AbsBead      point at parameter t of the named curve.
  BSubCurve    the part of the beads' parent curve between the two bead parameters.
  Line         straight segment between two points.
  RevSurf      profile revolved about the axis line (first axis point to second) through the
               stated angle range, right-handed about the axis direction. The range 270..360 deg
               turns the profile (in the plane y = 1, x < 0) through the quadrant x <= 0, y >= 1,
               which is what the deck's Extents line (y = +-2.513) requires.
  EdgeSnake    edge n of a RevSurf; edges run: 1 = profile at the start angle, 2 = revolved top
               end, 3 = profile at the end angle, 4 = revolved bottom end. Edge 3 is the
               profile in the plane y = 1.
  ProjCurve    orthogonal projection of a curve onto the plane Y=0.
  RuledSurf    straight lines joining equal parameters of the two curves.
  Symmetry     'x y' adds the images mirrored in x = 0, in y = 0 and in both.

Usage: python3 tests/reference/c1_reference.py [--deck PATH] z1 [z2 ...]
Prints a JSON list of {"z", "area", "Ixx", "Iyy", "x_half_width"}.
"""
import json
import sys

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import brentq


def parse_deck(path):
    text = open(path).read()
    symmetry = []
    for line in text.splitlines():
        if line.startswith("Symmetry:"):
            symmetry = line.split(":", 1)[1].split()
    body = text.split("BeginModel;", 1)[1].split("EndModel;", 1)[0]
    entities = {}
    for stmt in body.split(";"):
        stmt = stmt.strip().replace("{", " { ").replace("}", " } ")
        if not stmt:
            continue
        tokens = stmt.split()
        slash = tokens.index("/")
        kind, name = tokens[0], tokens[1]
        entities[name] = (kind, tokens[slash + 1:])
    return symmetry, entities


class Model:
    def __init__(self, entities):
        self.ent = entities
        self.cache = {}

    def point(self, name):
        if name == "*":
            return np.zeros(3)
        key = ("pt", name)
        if key not in self.cache:
            self.cache[key] = self._point(name)
        return self.cache[key]

    def _point(self, name):
        kind, a = self.ent[name]
        if kind == "FramePoint":
            return self.point(a[1]) + np.array([float(v) for v in a[3:6]])
        if kind == "MirrPoint":
            p = self.point(a[0]).copy()
            plane = a[1].lstrip("*")
            axis = {"X": 0, "Y": 1, "Z": 2}[plane[0]]
            level = float(plane.split("=")[1])
            p[axis] = 2 * level - p[axis]
            return p
        if kind == "AbsBead":
            return self.curve(a[0])(float(a[1]))
        raise ValueError("not a point entity: %s %s" % (kind, name))

    def curve(self, name):
        key = ("cv", name)
        if key not in self.cache:
            self.cache[key] = self._curve(name)
        return self.cache[key]

    def _curve(self, name):
        kind, a = self.ent[name]
        if kind == "BCurve":
            degree = int(a[1])
            pts = np.array([self.point(n) for n in a[a.index("{") + 1:a.index("}")]])
            n = len(pts)
            interior = [k / (n - degree) for k in range(1, n - degree)]
            knots = np.r_[np.zeros(degree + 1), interior, np.ones(degree + 1)]
            spline = BSpline(knots, pts, degree)
            return lambda t: spline(min(max(t, 0.0), 1.0))
        if kind == "Arc":
            p_a, centre, p_b = (self.point(n) for n in a[2:5])
            ua = p_a - centre
            radius = np.linalg.norm(ua)
            ua = ua / radius
            vb = p_b - centre
            vb = vb - np.dot(vb, ua) * ua
            wb = vb / np.linalg.norm(vb)
            sweep = np.arccos(np.clip(np.dot(p_b - centre, ua) / radius, -1, 1))
            return lambda t: centre + radius * (np.cos(t * sweep) * ua + np.sin(t * sweep) * wb)
        if kind == "PolyCurve2":
            parts = [self.curve(n) for n in a[a.index("{") + 1:a.index("}")]]
            m = len(parts)

            def poly(t):
                k = min(int(t * m), m - 1)
                return parts[k](t * m - k)
            return poly
        if kind == "BSubCurve":
            bead_a, bead_b = a[a.index("{") + 1:a.index("}")]
            parent = self.ent[bead_a][1][0]
            assert self.ent[bead_b][1][0] == parent, "beads lie on different curves"
            t0, t1 = float(self.ent[bead_a][1][1]), float(self.ent[bead_b][1][1])
            base = self.curve(parent)
            return lambda s: base(t0 + s * (t1 - t0))
        if kind == "Line":
            p0, p1 = self.point(a[1]), self.point(a[2])
            return lambda t: p0 + t * (p1 - p0)
        if kind == "EdgeSnake":
            surf = self.ent[a[2]]
            assert surf[0] == "RevSurf" and int(a[1]) == 3
            profile = self.curve(surf[1][1])
            return self._revolve(profile, surf[1][2], float(surf[1][4]))
        if kind == "ProjCurve":
            base = self.curve(a[1])
            axis = {"X": 0, "Y": 1, "Z": 2}[a[2].lstrip("*")[0]]

            def proj(t):
                p = base(t).copy()
                p[axis] = float(a[2].split("=")[1])
                return p
            return proj
        raise ValueError("not a curve entity: %s %s" % (kind, name))

    def _revolve(self, profile, axis_name, angle_deg):
        axis = self.ent[axis_name]
        a0, a1 = self.point(axis[1][1]), self.point(axis[1][2])
        k = (a1 - a0) / np.linalg.norm(a1 - a0)
        ang = np.deg2rad(angle_deg)

        def turned(t):
            v = profile(t) - a0
            return a0 + v * np.cos(ang) + np.cross(k, v) * np.sin(ang) + k * np.dot(k, v) * (1 - np.cos(ang))
        return turned


def z_roots(curve, z, n=4000):
    ts = np.linspace(0.0, 1.0, n + 1)
    zs = np.array([curve(t)[2] - z for t in ts])
    roots = []
    for i in range(n):
        if zs[i] == 0.0:
            roots.append(ts[i])
        elif zs[i] * zs[i + 1] < 0:
            roots.append(brentq(lambda t: curve(t)[2] - z, ts[i], ts[i + 1], xtol=1e-15, rtol=8.9e-16))
    if zs[n] == 0.0:
        roots.append(1.0)
    return roots


class Piece:
    """A boundary piece s in [0, 1] -> (x, y), with its analytic derivative d/ds."""

    def __init__(self, point, deriv):
        self.point = point
        self.deriv = deriv

    def flipped(self, fx, fy):
        f = np.array([fx, fy])
        return Piece(lambda s: self.point(s) * f, lambda s: self.deriv(s) * f)


def section_pieces(model, symmetry, z):
    pieces = []
    for name, (kind, a) in model.ent.items():
        if kind == "RevSurf":
            profile = model.curve(a[1])
            axis = model.ent[a[2]]
            a0, a1 = model.point(axis[1][1]), model.point(axis[1][2])
            k = (a1 - a0) / np.linalg.norm(a1 - a0)
            lo, hi = np.deg2rad(float(a[3])), np.deg2rad(float(a[4]))
            for t in z_roots(profile, z):
                v = profile(t) - a0
                kxv, kdv = np.cross(k, v), np.dot(k, v)

                def point(s, v=v, kxv=kxv, kdv=kdv, a0=a0, k=k, lo=lo, hi=hi):
                    ang = lo + s * (hi - lo)
                    return (a0 + v * np.cos(ang) + kxv * np.sin(ang) + k * kdv * (1 - np.cos(ang)))[:2]

                def deriv(s, v=v, kxv=kxv, kdv=kdv, k=k, lo=lo, hi=hi):
                    ang = lo + s * (hi - lo)
                    return ((hi - lo) * (-v * np.sin(ang) + kxv * np.cos(ang) + k * kdv * np.sin(ang)))[:2]
                pieces.append(Piece(point, deriv))
        elif kind == "RuledSurf":
            c1, c2 = model.curve(a[1]), model.curve(a[2])
            for t in z_roots(c1, z):
                p, q = c1(t), c2(t)
                pieces.append(Piece(lambda s, p=p, q=q: (p + s * (q - p))[:2],
                                    lambda s, p=p, q=q: (q - p)[:2]))
    flips = [(1.0, 1.0)]
    if "x" in symmetry:
        flips += [(-fx, fy) for fx, fy in flips]
    if "y" in symmetry:
        flips += [(fx, -fy) for fx, fy in flips]
    return [piece.flipped(fx, fy) for piece in pieces for fx, fy in flips]


def chain(pieces):
    """Order and orient the pieces into one closed loop by matching end points."""
    ends = [(p.point(0.0), p.point(1.0)) for p in pieces]
    used = [False] * len(pieces)
    order = [(0, False)]
    used[0] = True
    current = ends[0][1]
    for _ in range(len(pieces) - 1):
        best = None
        for i in range(len(pieces)):
            if used[i]:
                continue
            for rev, e in ((False, ends[i][0]), (True, ends[i][1])):
                d = np.linalg.norm(e - current)
                if best is None or d < best[0]:
                    best = (d, i, rev)
        d, i, rev = best
        if d > 1e-9:
            raise RuntimeError("section pieces do not join: gap %g" % d)
        used[i] = True
        order.append((i, rev))
        current = ends[i][0] if rev else ends[i][1]
    if np.linalg.norm(current - ends[0][0]) > 1e-9:
        raise RuntimeError("section loop does not close")
    return [(pieces[i], rev) for i, rev in order]


def loop_properties(loop):
    gx, gw = np.polynomial.legendre.leggauss(32)
    s = 0.5 * (gx + 1.0)
    w = 0.5 * gw
    area = ixx = iyy = 0.0
    xmax = 0.0
    for piece, rev in loop:
        ss = 1.0 - s if rev else s
        pts = np.array([piece.point(v) for v in ss])
        d = np.array([piece.deriv(v) for v in ss])
        if rev:
            d = -d
        x, y, dy = pts[:, 0], pts[:, 1], d[:, 1]
        # Green: A = closed integral of x dy; Ixx = integral of x y^2 dy; Iyy = integral of x^3/3 dy.
        area += np.sum(w * x * dy)
        ixx += np.sum(w * x * y * y * dy)
        iyy += np.sum(w * x ** 3 / 3.0 * dy)
        xmax = max(xmax, np.max(np.abs(x)))
    sign = 1.0 if area >= 0 else -1.0
    return sign * area, sign * ixx, sign * iyy, xmax


def main(argv):
    deck = "Input/C1.ms2"
    zs = []
    i = 1
    while i < len(argv):
        if argv[i] == "--deck":
            deck = argv[i + 1]
            i += 2
        else:
            zs.append(float(argv[i]))
            i += 1
    symmetry, entities = parse_deck(deck)
    model = Model(entities)
    result = []
    for z in zs:
        loop = chain(section_pieces(model, symmetry, z))
        area, ixx, iyy, xmax = loop_properties(loop)
        result.append({"z": z, "area": area, "Ixx": ixx, "Iyy": iyy, "x_half_width": xmax})
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main(sys.argv)
