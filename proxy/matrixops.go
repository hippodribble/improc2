package proxy

import (
	"log"

	"gonum.org/v1/gonum/mat"
)

type ImageMatrix struct {
	*mat.Dense
}

func NewImageMatrix(m *mat.Dense) ImageMatrix {
	return ImageMatrix{Dense: m}
}

// shifts an image by an integer number of pixels in x and y
//
//	dc,dr - shifts to apply - horizontal then vertical shifts
func (imat ImageMatrix) TranslatePixel(dc, dr int) ImageMatrix {
	// f := *im
	r, c := imat.Dims()

	var I, J, i, j int

	defer func() {
		if reco := recover(); reco != nil {
			log.Println("Recovered in f", i, j, I, J, r, c)
		}
	}()

	shifted := mat.NewDense(r, c, nil)
	for j = 0; j < r; j++ {
		J = j - dr

		for i = 0; i < c; i++ {
			I = i - dc
			if I > -1 && J > -1 && I < c && J < r {
				shifted.Set(j, i, imat.At(J, I))
			}
		}
	}
	return NewImageMatrix(shifted)
}

func (im ImageMatrix) InterpolateBilinear(row, col float64) float64 {
	row0 := int(row)
	row1 := row0 + 1
	col0 := int(col)
	col1 := col0 + 1
	dRow := row - float64(row0)
	dCol := col - float64(col0)
	defer func() {
		if r := recover(); r != nil {
			log.Println(row, col, row0, row1, col0, col1, dRow, dCol)
		}
	}()
	a := im.At(row0, col0)
	b := im.At(row1, col0)
	c := im.At(row0, col1)
	d := im.At(row1, col1)
	E := a + (b-a)*dRow
	F := c + (d-c)*dRow
	G := E + (F-E)*dCol
	return G
}

// TODO - finish this and create the bilinear interpolator
func (imat ImageMatrix) TranslateSubpixelSpaceDomainBilinear(dc, dr float64) ImageMatrix {

	r, c := imat.Dims()
	g:=ImageMatrix{Dense: mat.NewDense(r, c, nil)}
	
	var I, J float64
	for j := 0; j < r; j++ {
		for i := 0; i < c; i++ {
			I = float64(i) - dc
			J = float64(j) - dr
			if I >= 0 && I <= float64(c)-2 && J >= 0 && J <= float64(r)-2 {
				g.Set(j, i, imat.InterpolateBilinear(J, I))
			}
		}
	}
	return g
}

// returns image of gradient across columns (horizontally) within a matrix
func (imat ImageMatrix) ColumnGradient() ImageMatrix {

	r, c := imat.Dims()
	d := mat.NewDense(r, c, nil)
	for j := 0; j < r; j++ {
		for i := 1; i < c-1; i++ {
			d.Set(j, i, (imat.At(j, i+1)-imat.At(j, i-1))/2)
		}
	}
	return NewImageMatrix(d)
}

// returns image of gradient across rows (vertically) within a matrix
func (imat ImageMatrix) RowGradient() ImageMatrix {

	r, c := imat.Dims()
	d := mat.NewDense(r, c, nil)
	for i := 0; i < c; i++ {
		for j := 1; j < r-1; j++ {
			d.Set(j, i, (imat.At(j+1, i)-imat.At(j-1, i))/2)
		}
	}
	return NewImageMatrix(d)
}

// Horizontal gradient at an interpolated point in a matrix
//
// This corresponds to a single row r, which can be fractional, ie between discrete rows
// the gradient is returned across the specified fractional column
//
// There is no bounds checking, as this should be ensured before calling. This
// helps to improve speed.
func (imat *ImageMatrix) ColumnGradientAtPoint(r, c float64) float64 {
	return (imat.InterpolateBilinear(r, c+1) - imat.InterpolateBilinear(r, c-1)) / 2
}

func (imat *ImageMatrix) RowGradientAtPoint(r, c float64) float64 {
	return (imat.InterpolateBilinear(r+1, c) - imat.InterpolateBilinear(r-1, c)) / 2
}


func RescaleMatrixTo256(matrix mat.Matrix) mat.Dense {

	r, c := matrix.Dims()
	m := mat.NewDense(r, c, nil)
	vmin := 1e+12
	vmax := -1e+12

	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if matrix.At(i, j) < vmin {
				vmin = matrix.At(i, j)
			}
			if matrix.At(i, j) > vmax {
				vmax = matrix.At(i, j)
			}
		}
	}

	vrange := vmax - vmin
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			m.Set(i, j, (matrix.At(i, j)-vmin)/vrange*255)
		}
	}

	R, C := m.Dims()
	log.Println("Rescale to 256", R, C, "vs", r, c)
	
	return *m
}