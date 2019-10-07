// port of https://github.com/fontello/cubic2quad
package cubic2quad

import (
	"math"
	"sort"
)

func CalcPowerCoefficients(p1, c1, c2, p2 Point) []Point {
	// point(t) = p1*(1-t)^3 + c1*t*(1-t)^2 + c2*t^2*(1-t) + p2*t^3 = a*t^3 + b*t^2 + c*t + d
	// for each t value, so
	// a = (p2 - p1) + 3 * (c1 - c2)
	// b = 3 * (p1 + c2) - 6 * c1
	// c = 3 * (c1 - p1)
	// d = p1
	a := p2.Sub(p1).Add(c1.Sub(c2).Mul(3.0))
	b := p1.Add(c2).Mul(3.0).Sub(c1.Mul(6.0))
	c := c1.Sub(p1).Mul(3.0)
	d := p1
	return []Point{a, b, c, d}
}

func CalcPoint(a, b, c, d Point, t float64) Point {
	// a*t^3 + b*t^2 + c*t + d = ((a*t + b)*t + c)*t + d
	return a.Mul(t).Add(b).Mul(t).Add(c).Mul(t).Add(d)
}

func CalcPointQuad(a, b, c Point, t float64) Point {
	// a*t^2 + b*t + c = (a*t + b)*t + c
	return a.Mul(t).Add(b).Mul(t).Add(c)
}

func CalcPointDerivative(a, b, c, d Point, t float64) Point {
	// d/dt[a*t^3 + b*t^2 + c*t + d] = 3*a*t^2 + 2*b*t + c = (3*a*t + 2*b)*t + c
	return a.Mul(3.0 * t).Add(b.Mul(2.0)).Mul(t).Add(c)
}

func QuadSolve(a, b, c float64) []float64 {
	// a*x^2 + b*x + c = 0
	if a == 0.0 {
		if b == 0.0 {
			return []float64{}
		}
		return []float64{-c / b}
	}
	var D = b*b - 4.0*a*c
	if D < 0.0 {
		return []float64{}
	} else if D == 0.0 {
		return []float64{-b / (2.0 * a)}
	}
	var DSqrt = math.Sqrt(D)
	return []float64{(-b - DSqrt) / (2.0 * a), (-b + DSqrt) / (2.0 * a)}
}

func CubicRoot(x float64) float64 {
	if x < 0.0 {
		return -math.Pow(-x, 1.0/3.0)
	}
	return math.Pow(x, 1.0/3.0)
}

func CubicSolve(a, b, c, d float64) []float64 {
	// a*x^3 + b*x^2 + c*x + d = 0
	if a == 0.0 {
		return QuadSolve(b, c, d)
	}
	// solve using Cardan's method, which is described in paper of R.W.D. Nickals
	// http://www.nickalls.org/dick/papers/maths/cubic1993.pdf (doi:10.2307/3619777)
	xn := -b / (3.0 * a)                        // point of symmetry x coordinate
	yn := ((a*xn+b)*xn+c)*xn + d                // point of symmetry y coordinate
	deltaSq := (b*b - 3.0*a*c) / (9.0 * a * a)  // delta^2
	hSq := 4.0 * a * a * math.Pow(deltaSq, 3.0) // h^2
	D3 := yn*yn - hSq
	if D3 > 0.0 { // 1 real root
		var D3Sqrt = math.Sqrt(D3)
		return []float64{xn + CubicRoot((-yn+D3Sqrt)/(2.0*a)) + CubicRoot((-yn-D3Sqrt)/(2.0*a))}
	} else if D3 == 0.0 { // 2 real roots
		delta1 := CubicRoot(yn / (2.0 * a))
		return []float64{xn - 2.0*delta1, xn + delta1}
	}
	// 3 real roots
	var theta = math.Acos(-yn/math.Sqrt(hSq)) / 3.0
	var delta = math.Sqrt(deltaSq)
	return []float64{
		xn + 2.0*delta*math.Cos(theta),
		xn + 2.0*delta*math.Cos(theta+math.Pi*2.0/3.0),
		xn + 2.0*delta*math.Cos(theta+math.Pi*4.0/3.0),
	}
}

func MinDistanceToQuad(point, p1, c1, p2 Point) float64 {
	// f(t) = (1-t)^2 * p1 + 2*t*(1 - t) * c1 + t^2 * p2 = a*t^2 + b*t + c, t in [0, 1],
	// a = p1 + p2 - 2 * c1
	// b = 2 * (c1 - p1)
	// c = p1; a, b, c are vectors because p1, c1, p2 are vectors too
	// The distance between given point and quadratic curve is equal to
	// sqrt((f(t) - point)^2), so these expression has zero derivative by t at points where
	// (f'(t), (f(t) - point)) = 0.
	// Substituting quadratic curve as f(t) one could obtain a cubic equation
	// e3*t^3 + e2*t^2 + e1*t + e0 = 0 with following coefficients:
	// e3 = 2 * a^2
	// e2 = 3 * a*b
	// e1 = (b^2 + 2 * a*(c - point))
	// e0 = (c - point)*b
	// One of the roots of the equation from [0, 1], or t = 0 or t = 1 is a value of t
	// at which the distance between given point and quadratic Bezier curve has minimum.
	// So to find the minimal distance one have to just pick the minimum value of
	// the distance on set {t = 0 | t = 1 | t is root of the equation from [0, 1] }.

	a := p1.Add(p2).Sub(c1.Mul(2.0))
	b := c1.Sub(p1).Mul(2.0)
	c := p1
	e3 := 2.0 * a.Sqr()
	e2 := 3.0 * a.Dot(b)
	e1 := (b.Sqr() + 2.0*a.Dot(c.Sub(point)))
	e0 := c.Sub(point).Dot(b)
	candidates := []float64{}
	for _, t := range CubicSolve(e3, e2, e1, e0) {
		if t > 0.0 && t < 1.0 {
			candidates = append(candidates, t)
		}
	}
	candidates = append(candidates, 0.0, 1.0)

	minDistance := 1e9
	for i := 0; i < len(candidates); i++ {
		distance := CalcPointQuad(a, b, c, candidates[i]).Sub(point).Dist()
		if distance < minDistance {
			minDistance = distance
		}
	}
	return minDistance
}

func ProcessSegment(a, b, c, d Point, t1, t2 float64) []Point {
	// Find a single control point for given segment of cubic Bezier curve
	// These control point is an interception of tangent lines to the boundary points
	// Let's denote that f(t) is a vector function of parameter t that defines the cubic Bezier curve,
	// f(t1) + f'(t1)*z1 is a parametric equation of tangent line to f(t1) with parameter z1
	// f(t2) + f'(t2)*z2 is the same for point f(t2) and the vector equation
	// f(t1) + f'(t1)*z1 = f(t2) + f'(t2)*z2 defines the values of parameters z1 and z2.
	// Defining fx(t) and fy(t) as the x and y components of vector function f(t) respectively
	// and solving the given system for z1 one could obtain that
	//
	//      -(fx(t2) - fx(t1))*fy'(t2) + (fy(t2) - fy(t1))*fx'(t2)
	// z1 = ------------------------------------------------------.
	//            -fx'(t1)*fy'(t2) + fx'(t2)*fy'(t1)
	//
	// Let's assign letter D to the denominator and note that if D = 0 it means that the curve actually
	// is a line. Substituting z1 to the equation of tangent line to the point f(t1), one could obtain that
	// cx = [fx'(t1)*(fy(t2)*fx'(t2) - fx(t2)*fy'(t2)) + fx'(t2)*(fx(t1)*fy'(t1) - fy(t1)*fx'(t1))]/D
	// cy = [fy'(t1)*(fy(t2)*fx'(t2) - fx(t2)*fy'(t2)) + fy'(t2)*(fx(t1)*fy'(t1) - fy(t1)*fx'(t1))]/D
	// where c = (cx, cy) is the control point of quadratic Bezier curve.

	f1 := CalcPoint(a, b, c, d, t1)
	f2 := CalcPoint(a, b, c, d, t2)
	f1_ := CalcPointDerivative(a, b, c, d, t1)
	f2_ := CalcPointDerivative(a, b, c, d, t2)

	D := -f1_.x*f2_.y + f2_.x*f1_.y
	if math.Abs(D) < 1e-8 {
		return []Point{f1, f1.Add(f2).Div(2.0), f2} // straight line segment
	}
	cx := (f1_.x*(f2.y*f2_.x-f2.x*f2_.y) + f2_.x*(f1.x*f1_.y-f1.y*f1_.x)) / D
	cy := (f1_.y*(f2.y*f2_.x-f2.x*f2_.y) + f2_.y*(f1.x*f1_.y-f1.y*f1_.x)) / D
	return []Point{f1, Point{x: cx, y: cy}, f2}
}

func IsSegmentApproximationClose(a, b, c, d Point, tmin, tmax float64, p1, c1, p2 Point, errorBound float64) bool {
	// a,b,c,d define cubic curve
	// tmin, tmax are boundary points on cubic curve
	// p1, c1, p2 define quadratic curve
	// errorBound is maximum allowed distance
	// Try to find maximum distance between one of N points segment of given cubic
	// and corresponding quadratic curve that estimates the cubic one, assuming
	// that the boundary points of cubic and quadratic points are equal.
	//
	// The distance calculation method comes from Hausdorff distance defenition
	// (https://en.wikipedia.org/wiki/Hausdorff_distance), but with following simplifications
	// * it looks for maximum distance only for finite number of points of cubic curve
	// * it doesn't perform reverse check that means selecting set of fixed points on
	//   the quadratic curve and looking for the closest points on the cubic curve
	// But this method allows easy estimation of approximation error, so it is enough
	// for practical purposes.

	n := 10 // number of points + 1
	dt := (tmax - tmin) / float64(n)
	for t := tmin + dt; t < tmax-dt; t += dt { // don't check distance on boundary points
		// because they should be the same
		point := CalcPoint(a, b, c, d, t)
		if MinDistanceToQuad(point, p1, c1, p2) > errorBound {
			return false
		}
	}
	return true
}

func isApproximationClose(a, b, c, d Point, quadCurves [][]Point, errorBound float64) bool {
	dt := 1.0 / float64(len(quadCurves))
	for i := 0; i < len(quadCurves); i++ {
		p1 := quadCurves[i][0]
		c1 := quadCurves[i][1]
		p2 := quadCurves[i][2]
		if !IsSegmentApproximationClose(a, b, c, d, float64(i)*dt, float64(i+1)*dt, p1, c1, p2, errorBound) {
			return false
		}
	}
	return true
}

func FromFlatArray(points []float64) [][]Point {
	result := [][]Point{}
	segmentsNumber := (len(points) - 2) / 4
	for i := 0; i < segmentsNumber; i++ {
		result = append(result,
			[]Point{
				Point{x: points[4*i], y: points[4*i+1]},
				Point{x: points[4*i+2], y: points[4*i+3]},
				Point{x: points[4*i+4], y: points[4*i+5]},
			},
		)
	}
	return result
}

func ToFlatArray(quadsList [][]Point) []float64 {
	result := []float64{}
	result = append(result, quadsList[0][0].x)
	result = append(result, quadsList[0][0].y)
	for i := 0; i < len(quadsList); i++ {
		result = append(result, quadsList[i][1].x)
		result = append(result, quadsList[i][1].y)
		result = append(result, quadsList[i][2].x)
		result = append(result, quadsList[i][2].y)
	}
	return result
}

func IsApproximationClose(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y float64, quads []float64, errorBound float64) bool {
	// TODO: rewrite it in C-style and remove isApproximationClose
	pc := CalcPowerCoefficients(
		Point{p1x, p1y},
		Point{c1x, c1y},
		Point{c2x, c2y},
		Point{p2x, p2y},
	)
	return isApproximationClose(pc[0], pc[1], pc[2], pc[3], FromFlatArray(quads), errorBound)
}

/*
 * Split cubic bézier curve into two cubic curves, see details here:
 * https://math.stackexchange.com/questions/877725
 */
func SubdivideCubic(x1, y1, x2, y2, x3, y3, x4, y4, t float64) [][]float64 {
	var u = 1.0 - t
	var v = t

	var bx = x1*u + x2*v
	var sx = x2*u + x3*v
	var fx = x3*u + x4*v
	var cx = bx*u + sx*v
	var ex = sx*u + fx*v
	var dx = cx*u + ex*v

	var by = y1*u + y2*v
	var sy = y2*u + y3*v
	var fy = y3*u + y4*v
	var cy = by*u + sy*v
	var ey = sy*u + fy*v
	var dy = cy*u + ey*v

	return [][]float64{
		[]float64{x1, y1, bx, by, cx, cy, dx, dy},
		[]float64{dx, dy, ex, ey, fx, fy, x4, y4},
	}
}

func ByNumber(x, y float64) float64 {
	return x - y
}

type sortByNumber []float64

func (s sortByNumber) Len() int {
	return len(s)
}

func (s sortByNumber) Less(i, j int) bool {
	return i < j
}

func (s sortByNumber) Swap(i, j int) {
	t := s[i]
	s[i] = s[j]
	s[j] = t
}

/*
 * Find inflection points on a cubic curve, algorithm is similar to this one:
 * http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html
 */
func SolveInflections(x1, y1, x2, y2, x3, y3, x4, y4 float64) []float64 {
	var p = -(x4 * (y1 - 2.0*y2 + y3)) + x3*(2.0*y1-3.0*y2+y4) + x1*(y2-2.0*y3+y4) - x2*(y1-3.0*y3+2.0*y4)
	var q = x4*(y1-y2) + 3.0*x3*(-y1+y2) + x2*(2.0*y1-3.0*y3+y4) - x1*(2.0*y2-3.0*y3+y4)
	var r = x3*(y1-y2) + x1*(y2-y3) + x2*(-y1+y3)

	result := []float64{}
	for _, t := range QuadSolve(p, q, r) {
		if t > 1e-8 && t < 1-1e-8 {
			result = append(result, t)
		}
	}
	sort.Sort(sortByNumber(result))
	return result
}

/*
 * Approximate cubic Bezier curve defined with base points p1, p2 and control points c1, c2 with
 * with a few quadratic Bezier curves.
 * The function uses tangent method to find quadratic approximation of cubic curve segment and
 * simplified Hausdorff distance to determine number of segments that is enough to make error small.
 * In general the method is the same as described here: https://fontforge.github.io/bezier.html.
 */
func cubicToQuad(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, errorBound float64) []float64 {
	var p1 = Point{p1x, p1y}
	var c1 = Point{c1x, c1y}
	var c2 = Point{c2x, c2y}
	var p2 = Point{p2x, p2y}
	var pc = CalcPowerCoefficients(p1, c1, c2, p2)
	var a = pc[0]
	var b = pc[1]
	var c = pc[2]
	var d = pc[3]

	var approximation [][]Point
	for segmentsCount := 1; segmentsCount <= 8; segmentsCount++ {
		approximation = [][]Point{}
		for t := 0.0; t < 1.0; t += 1.0 / float64(segmentsCount) {
			approximation = append(approximation, ProcessSegment(a, b, c, d, t, t+1.0/float64(segmentsCount)))
		}
		if segmentsCount == 1 && (approximation[0][1].Sub(p1).Dot(c1.Sub(p1)) < 0.0 ||
			approximation[0][1].Sub(p2).Dot(c2.Sub(p2)) < 0.0) {
			// approximation concave, while the curve is convex (or vice versa)
			continue
		}
		if isApproximationClose(a, b, c, d, approximation, errorBound) {
			break
		}
	}
	return ToFlatArray(approximation)
}

/*
 * If this curve has any inflection points, split the curve and call
 * cubicToQuad function on each resulting curve.
 */
func CubicToQuad(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, errorBound float64) []float64 {
	var inflections = SolveInflections(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y)

	if len(inflections) == 0 {
		return cubicToQuad(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, errorBound)
	}

	result := []float64{}
	var curve = []float64{p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y}
	var prevPoint = 0.0

	for inflectionIdx := 0; inflectionIdx < len(inflections); inflectionIdx++ {
		split := SubdivideCubic(
			curve[0], curve[1], curve[2], curve[3],
			curve[4], curve[5], curve[6], curve[7],
			// we make a new curve, so adjust inflection point accordingly
			1.0-(1.0-inflections[inflectionIdx])/(1.0-prevPoint),
		)

		quad := cubicToQuad(
			split[0][0], split[0][1], split[0][2], split[0][3],
			split[0][4], split[0][5], split[0][6], split[0][7],
			errorBound,
		)

		result = append(result, quad[0:len(quad)-2]...)
		curve = split[1]
		prevPoint = inflections[inflectionIdx]
	}

	quad := cubicToQuad(
		curve[0], curve[1], curve[2], curve[3],
		curve[4], curve[5], curve[6], curve[7],
		errorBound,
	)

	result = append(result, quad...)
	return result
}
