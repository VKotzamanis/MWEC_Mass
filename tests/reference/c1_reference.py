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

Outer profile (the section of surface1 in its meridian plane y = 1, curve6, and equally the edge of
the ruled flat sides): points, unit tangents and outward unit normals, from analytic curve
derivatives. The normal of the revolved patch at another angle is the profile normal turned about
the axis; it is also computed from S_s x S_theta and the two must agree. The flat sides surface2 are
planes containing the y direction, so their normal is the profile normal as well.

Normal offset by t (the wall): the profile moved inward along its outward normal. A point of the
raw offset curve belongs to the offset surface only if its distance to the outer walls (the profile and its mirror image in x = 0) equals t;
points inside the fold or closer to the opposite wall (distance < t) are discarded. Revolution and the ruling leave this planar
construction unchanged, so the half-width of the offset section at a height is the |x| of the valid
offset-curve crossings of that height.

Usage:
  python3 tests/reference/c1_reference.py [--deck PATH] z1 [z2 ...]
      JSON list of {"z", "area", "Ixx", "Iyy", "x_half_width"}.
  python3 tests/reference/c1_reference.py [--deck PATH] --profile N
      JSON {"profile": [N points {"s", "x", "y", "z", "tx", "tz", "nx", "nz",
      "revolved_normal_gap"}]}: s is the curve6 parameter, (tx, tz) the unit tangent, (nx, nz) the
      outward unit normal in the profile plane, revolved_normal_gap the length of the difference
      between the turned normal and S_s x S_theta (unit) at the middle of the angle range (null
      for points on the axis).
  python3 tests/reference/c1_reference.py [--deck PATH] --offset T z1 [z2 ...]
      JSON list of {"z", "t", "x_half_width_valid": [...], "x_half_width_raw": [...],
      "x_half_width": largest valid value or null}: crossings of the height z by the offset curve,
      all of them and those that survive the fold test.
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

    def deriv(self, name):
        """Analytic d/dt of a curve entity, t in [0, 1]."""
        key = ("dcv", name)
        if key not in self.cache:
            self.cache[key] = self._deriv(name)
        return self.cache[key]

    def _deriv(self, name):
        kind, a = self.ent[name]
        if kind == "BCurve":
            degree = int(a[1])
            pts = np.array([self.point(n) for n in a[a.index("{") + 1:a.index("}")]])
            n = len(pts)
            interior = [k / (n - degree) for k in range(1, n - degree)]
            knots = np.r_[np.zeros(degree + 1), interior, np.ones(degree + 1)]
            d = BSpline(knots, pts, degree).derivative()
            return lambda t: d(min(max(t, 0.0), 1.0))
        if kind == "Arc":
            p_a, centre, p_b = (self.point(n) for n in a[2:5])
            ua = p_a - centre
            radius = np.linalg.norm(ua)
            ua = ua / radius
            vb = p_b - centre
            vb = vb - np.dot(vb, ua) * ua
            wb = vb / np.linalg.norm(vb)
            sweep = np.arccos(np.clip(np.dot(p_b - centre, ua) / radius, -1, 1))
            return lambda t: radius * sweep * (-np.sin(t * sweep) * ua + np.cos(t * sweep) * wb)
        if kind == "PolyCurve2":
            parts = [self.deriv(n) for n in a[a.index("{") + 1:a.index("}")]]
            m = len(parts)

            def dpoly(t):
                k = min(int(t * m), m - 1)
                return m * parts[k](t * m - k)
            return dpoly
        if kind == "BSubCurve":
            bead_a, bead_b = a[a.index("{") + 1:a.index("}")]
            parent = self.ent[bead_a][1][0]
            t0, t1 = float(self.ent[bead_a][1][1]), float(self.ent[bead_b][1][1])
            base = self.deriv(parent)
            return lambda s: (t1 - t0) * base(t0 + s * (t1 - t0))
        if kind == "Line":
            p0, p1 = self.point(a[1]), self.point(a[2])
            return lambda t: p1 - p0
        raise ValueError("no derivative for entity: %s %s" % (kind, name))

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


def profile_name(entities):
    surfaces = [(n, a) for n, (k, a) in entities.items() if k == "RevSurf"]
    assert len(surfaces) == 1, "one RevSurf expected"
    return surfaces[0][0], surfaces[0][1][1], surfaces[0][1]


def profile_frame(model, curve_name, s):
    """Point, unit tangent and outward unit normal (x-z plane) of the profile at parameter s.

    The profile runs from the top on the axis, towards -x and down to the keel: counter-clockwise
    in the x-z plane, so the outward normal is the tangent turned clockwise, (tz, -tx).
    """
    p = model.curve(curve_name)(s)
    d = model.deriv(curve_name)(s)
    t = np.array([d[0], d[2]]) / np.hypot(d[0], d[2])
    return p, t, np.array([t[1], -t[0]])


def revolved_normal_gap(model, surf_args, curve_name, s):
    """|turned profile normal - unit(S_s x S_theta)| at the middle of the angle range, or None on the axis."""
    axis = model.ent[surf_args[2]]
    a0, a1 = model.point(axis[1][1]), model.point(axis[1][2])
    k = (a1 - a0) / np.linalg.norm(a1 - a0)
    lo, hi = np.deg2rad(float(surf_args[3])), np.deg2rad(float(surf_args[4]))
    mid = 0.5 * (lo + hi)
    p, t, n = profile_frame(model, curve_name, s)
    d = model.deriv(curve_name)(s)
    v = p - a0

    def turn(w, ang):
        return w * np.cos(ang) + np.cross(k, w) * np.sin(ang) + k * np.dot(k, w) * (1 - np.cos(ang))
    # Profile plane: the profile sits at the end angle hi.
    n3 = np.array([n[0], 0.0, n[1]])
    turned = turn(n3, mid - hi)
    ss = turn(d, mid)
    stheta = np.cross(k, turn(v, mid))
    cross = np.cross(ss, stheta)
    if np.linalg.norm(stheta) < 1e-9:
        return None  # the profile point lies on the axis: S_theta vanishes, no normal from S_s x S_theta
    cross = cross / np.linalg.norm(cross)
    # S_s x S_theta may point either way; compare with the sign that agrees with the turned normal.
    if np.dot(cross, turned) < 0:
        cross = -cross
    return float(np.linalg.norm(turned - cross))


def point_segment_distance(points, poly):
    a, b = poly[:-1], poly[1:]
    ab = b - a
    denom = np.sum(ab * ab, axis=1)
    out = np.empty(len(points))
    for i, q in enumerate(points):
        u = np.clip(np.sum((q - a) * ab, axis=1) / denom, 0.0, 1.0)
        c = a + u[:, None] * ab
        out[i] = np.sqrt(np.min(np.sum((q - c) ** 2, axis=1)))
    return out


def offset_half_widths(model, curve_name, t_wall, zs, n_dense=20001, n_scan=4001):
    """Crossings of each height z by the inward normal offset of the profile, with the fold test."""
    dense_s = np.linspace(0.0, 1.0, n_dense)
    half = np.array([[model.curve(curve_name)(s)[0], model.curve(curve_name)(s)[2]] for s in dense_s])
    # The deck's symmetry 'x' puts the opposite wall at the mirror image, so a point of the void
    # must also be at least t from that wall (the two walls of a slender section never meet).
    dense = np.vstack([half, half[::-1] * np.array([-1.0, 1.0])])

    def offset_point(s):
        p, _, n = profile_frame(model, curve_name, s)
        return np.array([p[0], p[2]]) - t_wall * n

    scan = np.linspace(0.0, 1.0, n_scan)
    zoff = np.array([offset_point(s)[1] for s in scan])
    result = []
    for z in zs:
        raw = []
        for i in range(n_scan - 1):
            f0, f1 = zoff[i] - z, zoff[i + 1] - z
            if f0 == 0.0 or f0 * f1 < 0:
                root = brentq(lambda s: offset_point(s)[1] - z, scan[i], scan[i + 1], xtol=1e-15, rtol=8.9e-16)
                raw.append(root)
        pts = np.array([offset_point(r) for r in raw]).reshape(-1, 2)
        dist = point_segment_distance(pts, dense) if len(raw) else np.array([])
        # The distance of a valid offset point to the outer walls is t to within the dense
        # polyline's chord error; points of the fold are closer than t by far more than that.
        steps = np.hypot(*np.diff(half, axis=0).T)
        chord = float(np.max(steps))
        valid = [abs(pts[k, 0]) for k in range(len(raw)) if dist[k] >= t_wall - chord]
        result.append({"z": z, "t": t_wall,
                       "x_half_width_valid": valid,
                       "x_half_width_raw": [abs(p[0]) for p in pts],
                       "x_half_width": max(valid) if valid else None})
    return result


def main(argv):
    deck = "Input/C1.ms2"
    zs = []
    n_profile = None
    t_wall = None
    i = 1
    while i < len(argv):
        if argv[i] == "--deck":
            deck = argv[i + 1]
            i += 2
        elif argv[i] == "--profile":
            n_profile = int(argv[i + 1])
            i += 2
        elif argv[i] == "--offset":
            t_wall = float(argv[i + 1])
            i += 2
        else:
            zs.append(float(argv[i]))
            i += 1
    symmetry, entities = parse_deck(deck)
    model = Model(entities)
    if n_profile is not None or t_wall is not None:
        _, curve_name, surf_args = profile_name(entities)
        if n_profile is not None:
            rows = []
            for s in np.linspace(0.0, 1.0, n_profile):
                p, t, n = profile_frame(model, curve_name, s)
                rows.append({"s": s, "x": p[0], "y": p[1], "z": p[2], "tx": t[0], "tz": t[1],
                             "nx": n[0], "nz": n[1],
                             "revolved_normal_gap": revolved_normal_gap(model, surf_args, curve_name, s)})
            print(json.dumps({"profile": rows}, indent=1))
        else:
            print(json.dumps(offset_half_widths(model, curve_name, t_wall, zs), indent=1))
        return
    result = []
    for z in zs:
        loop = chain(section_pieces(model, symmetry, z))
        area, ixx, iyy, xmax = loop_properties(loop)
        result.append({"z": z, "area": area, "Ixx": ixx, "Iyy": iyy, "x_half_width": xmax})
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main(sys.argv)
