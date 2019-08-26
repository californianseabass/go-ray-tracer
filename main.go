package main

import (
	"image"
	"image/color"
	"image/png"
	"os"
)

func main() {
	nx := 200
	ny := 100
	image := image.NewRGBA(image.Rect(0, 0, nx, ny))

	for j := ny - 1; j >= 0; j-- {
		for i := 0; i < nx; i++ {
			r := uint8(float64(i) / float64(nx) * float64(255))
			g := uint8(float64(j) / float64(ny) * float64(255))
			b := uint8(float64(0.2) * float64(255))
			color := color.RGBA{r, g, b, 255}
			image.SetRGBA(i, j, color)
		}
	}

	outfile, err := os.Create("test.png")
	if err != nil {
		panic(err)
	}
	defer outfile.Close()
	err = png.Encode(outfile, image)
	if err != nil {
		panic(err)
	}
}
