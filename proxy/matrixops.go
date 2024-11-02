package proxy

import (
	"errors"
	"log"
	"math"

	"github.com/mjibson/go-dsp/fft"
	"gonum.org/v1/gonum/mat"
)

type ImageMatrix struct {
	*mat.Dense
}

func NewImageMatrix(m *mat.Dense) ImageMatrix {
	return ImageMatrix{Dense: m}
}

func NewPSF(radius int, sigma float64) ImageMatrix {

	psf := mat.NewDense(2*radius+1, 2*radius+1, nil)
	sum := 0.0
	for i := -radius; i <= radius; i++ {
		for j := -radius; j <= radius; j++ {
			a := math.Exp(-(math.Abs(float64(i)) + math.Abs(float64(j))) / 2 / sigma / sigma)
			psf.Set(i+radius, j+radius, a)
			sum += a
		}
	}
	divisor := 2*radius + 1
	divisor *= divisor

	for i := -radius; i <= radius; i++ {
		for j := -radius; j <= radius; j++ {
			psf.Set(i+radius, j+radius, psf.At(i+radius, j+radius)/float64(divisor))
		}
	}
	return NewImageMatrix(psf)

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
	g := ImageMatrix{Dense: mat.NewDense(r, c, nil)}

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
log.Printf("For picture, data range is %.1f - %.1f\n",vmin,vmax)
	vrange := vmax - vmin
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			m.Set(i, j, (matrix.At(i, j)-vmin)/vrange*255)
		}
	}

	// log.Println("Rescale to 256", R, C, "vs", r, c)

	return *m
}

func (f *ImageMatrix) Convolve(B ImageMatrix) ImageMatrix {

	var v float64

	R, C := f.Dims()
	r, _ := B.Dims()
	radius := r / 2
	g := mat.NewDense(R, C, nil)

	for i := 0; i < R; i++ {
		for j := 0; j < C; j++ {
			v = 0
			for k := -radius; k < radius+1; k++ {
				for l := -radius; l < radius+1; l++ {
					if i+k >= 0 && i+k < R && j+l >= 0 && j+l < C {
						v += f.At(i+k, j+l)
					}
				}
			}
			g.Set(i, j, v)
		}
	}
	return NewImageMatrix(g)
}

func (A *ImageMatrix) ToFloat2D() Float2D {

	r, _ := A.Dims()

	aA := make([][]float64, r)
	for i := 0; i < r; i++ {
		aA[i] = A.RawRowView(i)
	}
	return aA
}

func (A ImageMatrix) ConvolveFD(B ImageMatrix) (ImageMatrix, error) {
	// convert dense matrix to Float2D
	// and take spectra
	var product ImageMatrix

	var spA Complex2D = fft.FFT2Real(A.ToFloat2D())
	var spB Complex2D = fft.FFT2Real(B.ToFloat2D())
	// shift spectra
	spA = spA.Shift()
	spB = spB.Shift()
	// enlarge smaller spectrum to size of larger
	// ra, ca := spA.Dims()
	// rb, cb := spB.Dims()
	// log.Println("A dims:", ra, ca, "B dims:", rb, cb)
	// multiply complex spectra
	prodspec, err := spA.MultiplyElements(spB)
	if err != nil {
		return product, errors.New("error cross-multiplying spectra")
	}
	p := prodspec.Shift().IFFT().AsReal().ToDenseMatrix()
	return NewImageMatrix(&p), nil
}

func (g ImageMatrix) CalcRatio(f ImageMatrix) (ImageMatrix, error) {
	eta := 1.0e-30
	var G ImageMatrix
	r, c := g.Dims()
	R, C := f.Dims()
	if r != R || c != C {
		return G, errors.New("matrix dimensions do not match")
	}
	G = NewImageMatrix(mat.NewDense(r, c, nil))
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			G.Set(i, j, g.At(i, j)/(eta+f.At(i, j)))
		}
	}
	return G, nil
}

func (m ImageMatrix) Flip() ImageMatrix {
	r, c := m.Dims()
	var M mat.Dense = *mat.NewDense(r, c, nil)
	M.Copy(m)
	hr := r / 2
	hc := c / 2

	var a, b float64

	for i := 0; i < hr; i++ {
		for j := 0; j < c; j++ {
			a = M.At(r-1-i, j)
			b = M.At(i, j)
			M.Set(i, j, a)
			M.Set(r-1-i, j, b)
		}
	}
	for i := 0; i < r; i++ {
		for j := 0; j < hc; j++ {
			a = M.At(i, c-1-j)
			b = M.At(i, j)
			M.Set(i, j, a)
			M.Set(i, c-1-j, b)
		}
	}

	return NewImageMatrix(&M)
}

func (m *ImageMatrix) Cofactor(p, q int) ImageMatrix {
	n, _ := m.Dims()
	cofactor := mat.NewDense(n-1, n-1, nil)

	rowIndex, colIndex := 0, 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != p && j != q {
				cofactor.Set(rowIndex, colIndex, m.At(i, j))
				colIndex++
				if colIndex == n-1 {
					colIndex = 0
					rowIndex++
				}
			}
		}
	}
	return NewImageMatrix(cofactor)

}

func (m *ImageMatrix) Adjoint() ImageMatrix {
	n, _ := m.Dims()
	adj := NewImageMatrix(mat.NewDense(n, n, nil))

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			cofactor := m.Cofactor(i, j)
			sign := 1.0
			if (i+j)%2 != 0 {
				sign = -1.0
			}
			adj.Set(j, i, sign*mat.Det(cofactor))
		}
	}
	return adj
}

func (m *ImageMatrix) Pad(T, B, L, R int) *ImageMatrix {
	rows, cols := m.Dims()
	h := rows + T + B
	w := cols + L + R
	M := NewImageMatrix(mat.NewDense(h, w, nil))
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			M.Set(i+T, j+L, m.At(i, j))
		}
	}
	return &M
}

func (m *ImageMatrix) PadReflect(n int) *ImageMatrix {

	rows, cols := m.Dims()
	M := m.Pad(n, n, n, n)

	for i := 0; i < n; i++ {

		for col := 0; col < cols+2*n; col++ {
			M.Set(i, col, M.At(n*2-1-i, col))
			M.Set(rows+i+n, col, M.At(rows+n-i-1, col))
		}

		for row := 0; row < rows+2*n; row++ {
			M.Set(row, i, M.At(row, n*2-1-i))
			M.Set(row, cols+n+i, M.At(row, cols+n-i-1))
		}

	}

	return M
}

func (m *ImageMatrix) Trim(T, B, L, R int) *ImageMatrix {
	r, c := m.Dims()
	w := c - L - R
	h := r - T - B
	M := NewImageMatrix(mat.NewDense(h, w, nil))
	for i := 0; i < h; i++ {
		for j := 0; j < w; j++ {
			M.Set(i, j, m.At(i+T, j+L))
		}
	}
	return &M
}

func (m *ImageMatrix) TrimAll(n int) *ImageMatrix {
	return m.Trim(n, n, n, n)
}
