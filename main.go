package main

import (
	"image"
	"image/color"
	"image/png"
	"math"
	"math/rand"
	"os"
)

// camera.go
type Camera struct {
	Origin          Vector
	LowerLeftCorner Vector
	Horizontal      Vector
	Vertical        Vector
}

func (camera *Camera) getRay(u float64, v float64) *Ray {
	var direction = camera.LowerLeftCorner.Add(camera.Horizontal.Scale(u).Add(camera.Vertical.Scale(v)).Minus(&camera.Origin))
	return &Ray{
		camera.Origin,
		*direction,
	}
}

var MainCamera Camera = Camera{
	LowerLeftCorner: Vector{-2.0, -1.0, -1.0},
	Horizontal:      Vector{4.0, 0.0, 0.0},
	Vertical:        Vector{0.0, 2.0, 0.0},
	Origin:          Vector{0.0, 0.0, 0.0},
}

// hitable.go
type HitRecord struct {
	T        float64
	P        *Vector
	Normal   *Vector
	Material Material
}

type Hitable interface {
	hit(r *Ray, t_min float64, t_max float64, rec *HitRecord) bool
}

// metal.go
type Metal struct {
	albedo Vector
}

func (metal Metal) Scatter(ray *Ray, rec *HitRecord, attenuation *Vector, scattered *Ray) bool {
	reflected := reflect(*UnitVector(&ray.Direction), *rec.Normal)
	*scattered = Ray{*rec.P, *reflected}
	*attenuation = metal.albedo
	return true
}

// lambertian.go
type Lambertian struct {
	albedo Vector
}

func (lambertian Lambertian) Scatter(ray *Ray, rec *HitRecord, attenuation *Vector, scattered *Ray) bool {
	target := rec.P.Add(rec.Normal).Add(RandomInUnitSphere())
	*scattered = Ray{*rec.P, *target.Minus(rec.P)}
	*attenuation = lambertian.albedo
	return true
}

// hitable_list.go
type HitableList []Hitable

func (hs HitableList) hit(r *Ray, t_min float64, t_max float64, rec *HitRecord) bool {
	temp := HitRecord{}
	hitAnything := false
	closestSoFar := t_max
	for _, h := range hs {
		if h.hit(r, t_min, closestSoFar, &temp) {
			hitAnything = true
			closestSoFar = temp.T
			*rec = temp
		}
	}
	return hitAnything
}

type Material interface {
	Scatter(ray *Ray, rec *HitRecord, attenuation *Vector, scattered *Ray) bool
}

// sphere.go

type Sphere struct {
	Center   Vector
	Radius   float64
	Material Material
}

func (s *Sphere) hit(ray *Ray, t_min float64, t_max float64, rec *HitRecord) bool {
	oc := ray.Origin.Minus(&s.Center)
	a := ray.Direction.Dot(&ray.Direction)
	b := oc.Dot(&ray.Direction)
	c := oc.Dot(oc) - s.Radius*s.Radius
	discriminant := b*b - a*c
	if discriminant > 0 {
		temp := (-b - math.Sqrt(discriminant)) / a
		if temp < t_max && temp > t_min {
			rec.T = temp
			rec.P = ray.point_at_parameter(temp)
			rec.Normal = rec.P.Minus(&s.Center).Divide(s.Radius)
			rec.Material = s.Material
			return true
		}
		temp = (-b + math.Sqrt(discriminant)) / a
		if temp < t_max && temp > t_min {
			rec.T = temp
			rec.P = ray.point_at_parameter(temp)
			rec.Normal = rec.P.Minus(&s.Center).Divide(s.Radius)
			rec.Material = s.Material
			return true
		}
	}
	return false
}

func RandomInUnitSphere() *Vector {
	var vec Vector
	lengthOneVector := &Vector{1.0, 1.0, 1.0}
	for i := 0; true; i++ {
		v := Vector{rand.Float64(), rand.Float64(), rand.Float64()}
		vec = *v.Scale(2.0).Minus(lengthOneVector)
		if vec.SquaredLength() < 1.0 {
			break
		}
	}
	return &vec
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

func (v *Vector) Scale(t float64) *Vector {
	return &Vector{v.X * t, v.Y * t, v.Z * t}
}

func (v *Vector) Multiply(u *Vector) *Vector {
	return &Vector{v.X * u.X, v.Y * u.Y, v.Z * u.Z}
}

func (v *Vector) Divide(t float64) *Vector {
	return &Vector{v.X / t, v.Y / t, v.Z / t}
}

func UnitVector(v *Vector) *Vector {
	length := v.Length()
	return v.Divide(length)
}

func (v *Vector) Length() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

func (v *Vector) Dot(u *Vector) float64 {
	return v.X*u.X + v.Y*u.Y + v.Z*u.Z
}

func (v *Vector) SquaredLength() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

func reflect(v Vector, n Vector) *Vector {
	return v.Minus(n.Scale(v.Dot(&n) * 2))
}

// main.go

func backgroundPixel(ray *Ray) *Vector {
	unit_direction := UnitVector(&ray.Direction)
	// scale the y contribution to be between 0 and 1
	t := 0.5 * (unit_direction.Y + 1.0)
	// lerp = blended_value = (1-t)*start_value + t*end_value
	startValue := &Vector{1.0, 1.0, 1.0}
	endValue := &Vector{0.5, 0.7, 1.0}
	pixel := startValue.Scale(1.0 - t).Add(endValue.Scale(t))
	return pixel
}

func VectorToRGBA(pixel *Vector) *color.RGBA {
	// the sqrt is a gamma 2 correction
	r := uint8(float64(255.99) * math.Sqrt(pixel.X))
	g := uint8(float64(255.99) * math.Sqrt(pixel.Y))
	b := uint8(float64(255.99) * math.Sqrt(pixel.Z))
	return &color.RGBA{r, g, b, 255}
}

func pixelValue(ray *Ray, hs HitableList, depth int) *Vector {
	rec := &HitRecord{}
	isHit := hs.hit(ray, 0.001, math.MaxFloat64, rec)
	if isHit {
		scattered := &Ray{}
		attenuation := &Vector{}
		if depth < 50 && rec.Material.Scatter(ray, rec, attenuation, scattered) {
			return attenuation.Multiply(pixelValue(scattered, hs, depth+1))
		} else {
			return &Vector{0.0, 0.0, 0.0}
		}
	} else {
		return backgroundPixel(ray)
	}
}

var circleCenter Vector = Vector{0.0, 0.0, -1.0}
var circle Sphere = Sphere{circleCenter, 0.5, Lambertian{Vector{0.8, 0.3, 0.3}}}
var groundCenter Vector = Vector{0.0, -100.5, -1}
var ground = Sphere{groundCenter, 100.0, Lambertian{Vector{0.8, 0.8, 0.0}}}
var circleACenter Vector = Vector{1.0, 0.0, -1.0}
var circleA Sphere = Sphere{circleACenter, 0.5, Metal{Vector{0.8, 0.6, 0.2}}}
var circleBCenter Vector = Vector{-1.0, 0.0, -1.0}
var circleB Sphere = Sphere{circleBCenter, 0.5, Metal{Vector{0.8, 0.8, 0.8}}}

func itemsInScene() HitableList {
	hs := make([]Hitable, 4)
	hs[0] = &circle
	hs[1] = &ground
	hs[2] = &circleA
	hs[3] = &circleB
	return hs
}

func main() {
	nx := 600
	ny := 300
	nSamples := 100
	image := image.NewRGBA(image.Rect(0, 0, nx, ny))

	world := itemsInScene()

	for j := ny - 1; j >= 0; j-- {
		for i := 0; i < nx; i++ {
			pixel := &Vector{0, 0, 0}
			// take multiple samples to anti alias and blend boundaries
			for s := 0; s < nSamples; s++ {
				u := (float64(i) + rand.Float64()) / float64(nx) // x coordinate
				v := (float64(j) + rand.Float64()) / float64(ny) // y coordinte
				ray := MainCamera.getRay(u, v)
				pixel = pixel.Add(pixelValue(ray, world, 0))
			}
			pixel = pixel.Divide(float64(nSamples))
			color := VectorToRGBA(pixel)
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
