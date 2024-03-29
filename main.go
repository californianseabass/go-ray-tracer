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
	LensRadius      float64
	u               Vector
	v               Vector
	w               Vector
}

// change this for camera angles that aren't level to the horizontal plane
var vup Vector = Vector{0, 1, 0}

func randomInUnitDisk() Vector {
	p := Vector{math.MaxFloat64, math.MaxFloat64, math.MaxFloat64}
	for true {
		if p.Dot(&p) < 1.0 {
			break
		}
		randomVector := Vector{rand.Float64(), rand.Float64(), 0}
		p = *randomVector.Scale(2.0).Minus(&Vector{1, 1, 0})
	}
	return p
}

func (camera *Camera) New(lookFrom Vector, lookAt Vector, vfovDegrees float64, aspect float64, aperture float64, focusDist float64) {
	camera.LensRadius = aperture / 2
	theta := vfovDegrees * math.Pi / 180
	halfHeight := math.Tan(theta / 2)
	halfWidth := aspect * halfHeight
	camera.Origin = lookFrom

	w := UnitVector(lookFrom.Minus(&lookAt))
	u := UnitVector(vup.Cross(w))
	v := w.Cross(u)
	camera.w = *w
	camera.u = *u
	camera.v = *v
	camera.LowerLeftCorner = *lookFrom.Minus(u.Scale(halfWidth * focusDist)).Minus(v.Scale(halfHeight * focusDist)).Minus(w.Scale(focusDist))
	camera.Horizontal = *u.Scale(2 * halfWidth * focusDist)
	camera.Vertical = *v.Scale(2 * halfHeight * focusDist)
}

func (camera Camera) getRay(s float64, t float64) *Ray {
	rd := randomInUnitDisk()
	rd = *rd.Scale(camera.LensRadius)
	offset := camera.u.Scale(rd.X).Add(camera.v.Scale(rd.Y))
	var direction = camera.LowerLeftCorner.Add(camera.Horizontal.Scale(s)).Add(camera.Vertical.Scale(t)).Minus(&camera.Origin).Minus(offset)
	return &Ray{
		*camera.Origin.Add(offset),
		*direction,
	}
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

// dielectric.go
func schlick(cosine float64, refIdx float64) float64 {
	r0 := (1 - refIdx) / (1 + refIdx)
	r0 = r0 * r0
	return r0 + (1-r0)*math.Pow((1-cosine), 5)
}

type Dielectric struct {
	refIdx float64
}

func (de Dielectric) Scatter(ray *Ray, rec *HitRecord, attenuation *Vector, scattered *Ray) bool {
	outwardNormal := Vector{}
	refracted := &Vector{}
	reflected := reflect(ray.Direction, *rec.Normal)
	var niOverNt float64
	var cosine float64
	*attenuation = Vector{1.0, 1.0, 1.0}
	var reflectProb float64
	if ray.Direction.Dot(rec.Normal) > 0 {
		outwardNormal = *rec.Normal.Negate()
		niOverNt = de.refIdx
		cosine = de.refIdx * ray.Direction.Dot(rec.Normal) / ray.Direction.Length()
	} else {
		outwardNormal = *rec.Normal
		niOverNt = 1.0 / de.refIdx
		cosine = -1.0 * ray.Direction.Dot(rec.Normal) / ray.Direction.Length()
	}
	if Refract(ray.Direction, outwardNormal, niOverNt, refracted) {
		reflectProb = schlick(cosine, de.refIdx)
	} else {
		*scattered = Ray{*rec.P, *reflected}
		reflectProb = 1.0
	}
	if rand.Float64() < reflectProb {
		*scattered = Ray{*rec.P, *reflected}
	} else {
		*scattered = Ray{*rec.P, *refracted}
	}

	return true
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

func (v *Vector) Negate() *Vector {
	return &Vector{-v.X, -v.Y, -v.Z}
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

func (v *Vector) Cross(u *Vector) *Vector {
	return &Vector{v.Y*u.Z - v.Z*u.Y, -(v.X*u.Z - v.Z*u.X), v.X*u.Y - v.Y*u.X}
}

func (v *Vector) SquaredLength() float64 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

func reflect(v Vector, n Vector) *Vector {
	return v.Minus(n.Scale(v.Dot(&n) * 2))
}

func Refract(v Vector, n Vector, niOverNt float64, refracted *Vector) bool {
	unit := UnitVector(&v)
	dt := unit.Dot(&n)
	discriminant := 1.0 - niOverNt*niOverNt*(1-dt*dt)
	if discriminant > 0.0 {
		*refracted = *unit.Minus(n.Scale(dt)).Scale(niOverNt).Minus(n.Scale(math.Sqrt(discriminant)))
		return true
	} else {
		return false
	}
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
			temp := attenuation.Multiply(pixelValue(scattered, hs, depth+1))
			return temp
		} else {
			return &Vector{0.0, 0.0, 0.0}
		}
	} else {
		return backgroundPixel(ray)
	}
}

func itemsInScene() HitableList {
	hs := []Hitable{}
	hs = append(hs, &Sphere{Vector{0, -1000, 0}, 1000, Lambertian{Vector{0.5, 0.5, 0.5}}})
	for a := -11; a < 11; a++ {
		for b := -11; b < 11; b++ {
			chooseMaterial := rand.Float64()
			center := Vector{float64(a) + 0.9*rand.Float64(), 0.2, float64(b) + 0.9*rand.Float64()}
			if center.Minus(&Vector{4, 0.2, 0}).Length() > 0.9 {
				if chooseMaterial < 0.8 { // diffuse
					hs = append(hs, &Sphere{center, 0.2, Lambertian{Vector{rand.Float64(), rand.Float64(), rand.Float64()}}})
				} else if chooseMaterial < 0.95 { // metal
					hs = append(hs, &Sphere{center, 0.2, Metal{Vector{0.5 * (1 + rand.Float64()), 0.5 * (1 + rand.Float64()), 0.5 * (1 + rand.Float64())}}})
				} else { // dielectic
					hs = append(hs, &Sphere{center, 0.2, Dielectric{1.5}})
				}
			}
		}
	}
	hs = append(hs, &Sphere{Vector{0, 1, 0}, 1.0, Dielectric{1.5}})
	hs = append(hs, &Sphere{Vector{-4, 1, -0}, 1.0, Lambertian{Vector{0.4, 0.2, 0.1}}})
	hs = append(hs, &Sphere{Vector{4, 1, 0}, 1.0, Metal{Vector{0.7, 0.6, 0.5}}})
	return hs
}

func main() {
	MainCamera := &Camera{}
	nx := 600
	ny := 300
	nSamples := 100
	image := image.NewRGBA(image.Rect(0, 0, nx, ny))

	lookFrom := Vector{13, 2, 2}
	lookAt := Vector{0, 0, 0}
	distToFocus := lookFrom.Minus(&lookAt).Length()
	MainCamera.New(lookFrom, lookAt, 20, float64(nx)/float64(ny), 0.1, distToFocus)

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
