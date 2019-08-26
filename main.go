package main

import (
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
)

// hitable.go
type hit_record struct {
	t      float64
	p      Vector
	normal Vector
}

type hitable interface {
	hit(r *Ray, t_min float64, t_max float64, rec *hit_record) bool
}

// sphere.go

type Sphere struct {
	Center Vector
	Radius float64
}

func (s *Sphere) hit(r *Ray) bool {
	oc := r.Origin.Minus(&s.Center)
	a := r.Direction.Dot(&r.Direction)
	b := oc.Dot(&r.Direction) * 2.0
	c := oc.Dot(oc) - s.Radius*s.Radius
	discriminant := b*b - 4*a*c
	return discriminant > 0
}

// ray.go
type Ray struct {
	Origin    Vector
	Direction Vector
}

func (r *Ray) point_at_parameter(t float64) *Vector {
	return &Vector{
		r.Origin.X + r.Direction.X*t,
		r.Origin.Y + r.Direction.Y*t,
		r.Origin.Z + r.Direction.Z*t,
	}
}

// vector.go
type Vector struct {
	X float64
	Y float64
	Z float64
}

func (v *Vector) Add(u *Vector) *Vector {
	return &Vector{v.X + u.X, v.Y + u.Y, v.Z + u.Z}
}

func (v *Vector) Minus(u *Vector) *Vector {
	return &Vector{v.X - u.X, v.Y - u.Y, v.Z - u.Z}
}

func (v *Vector) Multiply(t float64) *Vector {
	return &Vector{v.X * t, v.Y * t, v.Z * t}
}

func (v *Vector) Divide(t float64) *Vector {
	return &Vector{v.X / t, v.Y / t, v.Z / t}
}

func UnitVector(v Vector) *Vector {
	length := v.Length()
	return v.Divide(length)
}

func (v *Vector) Length() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

func (v *Vector) Dot(u *Vector) float64 {
	return v.X*u.X + v.Y*u.Y + v.Z*u.Z
}

// main.go
var redCircleCenter Vector = Vector{0.0, 0.0, 1.0}
var redCircle *Sphere = &Sphere{redCircleCenter, 0.5}

func backgroundPixel(ray *Ray) *Vector {
	unit_direction := UnitVector(ray.Direction)
	// scale the y contribution to be between 0 and 1
	t := 0.5 * (unit_direction.Y + 1.0)
	// lerp = blended_value = (1-t)*start_value + t*end_value
	startValue := &Vector{1.0, 1.0, 1.0}
	endValue := &Vector{0.5, 0.7, 1.0}
	pixel := startValue.Multiply(1.0 - t).Add(endValue.Multiply(t))
	return pixel
}

func pixelValue(ray *Ray) *color.RGBA {
	var pixel *Vector
	if redCircle.hit(ray) {
		pixel = &Vector{1.0, 0.0, 0.0}
	} else {
		pixel = backgroundPixel(ray)
	}
	r := uint8(float64(255.99) * pixel.X)
	g := uint8(float64(255.99) * pixel.Y)
	b := uint8(float64(255.99) * pixel.Z)
	color := &color.RGBA{r, g, b, 255}
	return color
}

func main() {
	nx := 200
	ny := 100
	image := image.NewRGBA(image.Rect(0, 0, nx, ny))

	lower_left_corner := Vector{-2.0, -1.0, -1.0}
	horizontal := Vector{4.0, 0.0, 0.0}
	vertical := Vector{0.0, 2.0, 0.0}
	origin := Vector{0.0, 0.0, 0.0}

	for j := ny - 1; j >= 0; j-- {
		for i := 0; i < nx; i++ {
			u := float64(i) / float64(nx) // x coordinate
			v := float64(j) / float64(ny) // y coordinte

			direction := lower_left_corner.Add(horizontal.Multiply(u).Add(vertical.Multiply(v)))
			ray := &Ray{origin, *direction}

			color := pixelValue(ray)
			image.SetRGBA(i, int(math.Floor(math.Abs(float64(j-ny)))), *color)
		}
	}

	outfile, err := os.Create("output.png")
	if err != nil {
		panic(err)
	}
	defer outfile.Close()
	err = png.Encode(outfile, image)
	if err != nil {
		panic(err)
	}
}
