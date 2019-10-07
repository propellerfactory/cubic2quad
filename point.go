package cubic2quad

// port of https://github.com/fontello/cubic2quad
import "math"

type Point struct {
	x, y float64
}

func (p Point) Add(op Point) Point {
	return Point{
		x: p.x + op.x,
		y: p.y + op.y,
	}
}

func (p Point) Sub(op Point) Point {
	return Point{
		x: p.x - op.x,
		y: p.y - op.y,
	}
}

func (p Point) Mul(v float64) Point {
	return Point{
		x: p.x * v,
		y: p.y * v,
	}
}

func (p Point) Div(v float64) Point {
	return Point{
		x: p.x / v,
		y: p.y / v,
	}
}

func (p Point) Dist() float64 {
	return math.Sqrt(p.x*p.x + p.y*p.y)
}

func (p Point) Sqr() float64 {
	return p.x*p.x + p.y*p.y
}

func (p Point) Dot(op Point) float64 {
	return p.x*op.x + p.y*op.y
}
