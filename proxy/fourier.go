package proxy

import (
	"errors"
	"math"
	"math/cmplx"

	"github.com/mjibson/go-dsp/fft"
)

// wrapper for a 2D complex array
type Complex2D [][]complex128

func (c Complex2D) Spectrum() Complex2D {
	return fft.FFT2(c)
}

func (c Complex2D) IFFT() Complex2D {
	return fft.IFFT2(c)
}

func (c Complex2D) AsAmplitude() Float2D {
	var out Float2D = make([][]float64, len(c))
	for i := 0; i < len(c); i++ {
		out[i] = make([]float64, len(c[0]))
		for j := 0; j < len(c[0]); j++ {
			out[i][j] = cmplx.Abs(c[i][j])
		}
	}
	return out
}

func (c Complex2D) AsPhase() Float2D {
	var out Float2D = make([][]float64, len(c))
	for i := 0; i < len(c); i++ {
		out[i] = make([]float64, len(c[0]))
		for j := 0; j < len(c[0]); j++ {
			out[i][j] = cmplx.Phase(c[i][j])
		}
	}
	return out
}

func (c Complex2D) AsReal() Float2D {
	var out Float2D = make([][]float64, len(c))
	for i := 0; i < len(c); i++ {
		out[i] = make([]float64, len(c[0]))
		for j := 0; j < len(c[0]); j++ {
			out[i][j] = real(c[i][j])
		}
	}
	return out
}

func (c Complex2D) AsImaginary() Float2D {
	var out Float2D = make([][]float64, len(c))
	for i := 0; i < len(c); i++ {
		out[i] = make([]float64, len(c[0]))
		for j := 0; j < len(c[0]); j++ {
			out[i][j] = imag(c[i][j])
		}
	}
	return out
}

// TranslateFD applies a spatial delay in 2D to a Spectrum2D of an image, by
//   - creating an appropriate 2D phase-only filter, representing delay
//   - applying the filter to the data
func (c Complex2D) TranslateFD(dx, dy float64) Complex2D {

	w := len(c[0])
	h := len(c)

	outspec := make([][]complex128, h)
	phaseshiftspec := make([][]complex128, h)

	var f, g float64
	for j := 0; j < h; j++ {
		phaseshiftspec[j] = make([]complex128, w)
		outspec[j] = make([]complex128, w)
		jmod := float64(j + h/2)
		jmod = math.Mod(jmod, float64(h))
		jmod /= float64(h)
		f = -2 * math.Pi * dy * (jmod - 0.5)
		for i := 0; i < w; i++ {
			imod := float64(i + w/2)
			imod = math.Mod(imod, float64(w))
			imod /= float64(w)
			g = -2 * math.Pi * dx * (imod - 0.5)
			phaseshiftspec[j][i] = cmplx.Exp(complex(0, f+g))
		}
	}

	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			outspec[j][i] = c[j][i] * phaseshiftspec[j][i]
		}
	}
	return outspec
}

// returns the integer sample location of the peak of a Complex2D
func (c Complex2D) MaxArgs() (int, int) {
	h := len(c)
	w := len(c[0])
	var a float64
	maxamp := 0.0
	var xmax, ymax int

	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			a = cmplx.Abs(c[j][i])

			if a > maxamp {
				maxamp = a
				xmax = i
				ymax = j
			}
		}
	}
	if xmax > w/2 {
		xmax = xmax - w
	}
	if ymax > h/2 {
		ymax = ymax - h
	}
	return xmax, ymax
}


// gets the subsample peak of a Complex2D
func (c Complex2D) MaxPeakFromGaussianFit() (float64, float64) {
	// spec := *specp
	H := len(c)
	W := len(c[0])
	X, Y := c.MaxArgs()

	// get neighbours and middle
	M := math.Log10(cmplx.Abs(c[(Y+H)%H][(X+W)%W]))
	L := math.Log10(cmplx.Abs(c[(Y+H)%H][(X-1+W)%W]))
	R := math.Log10(cmplx.Abs(c[(Y+H)%H][(X+1+W)%W]))
	T := math.Log10(cmplx.Abs(c[(Y-1+H)%H][(X+W)%W]))
	B := math.Log10(cmplx.Abs(c[(Y+1+H)%H][(X+W)%W]))

	gaussX := (L - R) / (L + R - 2*M) / 2
	gaussY := (T - B) / (T + B - 2*M) / 2

	// fmt.Printf("GAUSSIAN PEAKS: %.2f , %.2f\n", gaussX, gaussY)

	// ShowCurrentOperation("")
	return gaussX + float64(X), gaussY + float64(Y)
}

// shifts the float image by half in x and y, to central DC in amplitude spectra etc.
func (c Complex2D) Shift() Complex2D {

	w := len(c[0])
	h := len(c)

	var cOut Complex2D = make([][]complex128, h)

	for j := 0; j < h; j++ {
		cOut[j] = make([]complex128, w)
	}

	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			J := (j + h/2) % h
			I := (i + w/2) % w
			cOut[J][I] = c[j][i]
		}
	}
	return cOut
}

func(c Complex2D)Dims()(rows,cols int){
	return len(c),len(c[0])
}

func(c Complex2D)MultiplyElements(a Complex2D)(Complex2D,error){
	R,C:=c.Dims()
	R2,C2:=a.Dims()
	if R!=R2 || C!=C2{
		return nil, errors.New("spectral dimensions do not match")
	}
	var out Complex2D=make([][]complex128,R)
	for i:=0;i<R;i++{
		out[i]=make([]complex128,C)
		for j:=0;j<C;j++{
			out[i][j]=c[i][j]*a[i][j]
		}
	}
	return out,nil
}