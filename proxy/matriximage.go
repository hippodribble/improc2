package proxy

import (
	"image"
	"image/color"

	"gonum.org/v1/gonum/mat"
)

type MatrixImage struct {
	matrix mat.Dense
}

func NewMatrixImage(matrix mat.Matrix) *MatrixImage {
	im := &MatrixImage{matrix: RescaleMatrixTo256(matrix)}
	return im
}

func (mi *MatrixImage) At(x, y int) color.Color {
	return color.Gray{uint8(mi.matrix.At(y, x))}
}

func (mi *MatrixImage) ColorModel() color.Model {
	return color.GrayModel
}

func (mi *MatrixImage) Bounds() image.Rectangle {
	r, c := mi.matrix.Dims()
	return image.Rect(0, 0, c, r)
}
